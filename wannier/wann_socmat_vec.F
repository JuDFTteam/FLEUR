!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter GrÃ¼nberg Institut, Forschungszentrum JÃ¼lich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_socmat_vec
      contains
      subroutine wann_socmat_vec(
     >               mpi,enpara,input,noco,atoms,
     >               lmaxd,natd,neigd,
     >               llod,jspd,
     >               jspins,irank,
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

      TYPE(t_mpi),INTENT(IN)      :: mpi

      TYPE(t_enpara),INTENT(IN)   :: enpara
      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_noco),INTENT(IN)     :: noco
      TYPE(t_atoms),INTENT(IN)    :: atoms

      integer, intent(in) :: lmaxd
      integer, intent(in) :: natd
      integer, intent(in) :: neigd

      integer, intent(in) :: llod
      integer, intent(in) :: jspd

      integer, intent(in) :: jspins
      integer, intent(in) :: irank

      real,    intent(in) :: vr(:,0:,:,:)    ! jmtd,0:nlhd,ntypd,jspd   
  
      complex, intent(in) :: acof(:,0:,:,:) !acof(noccbd,0:lmd,natd,jspd)
      complex, intent(in) :: bcof(:,0:,:,:) !bcof(noccbd,0:lmd,natd,jspd)
      complex, intent(in) :: chelp(-llod:,:,:,:,:) !chelp(-llod:llod,neigd,nlod,natd,jspd)

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

      nwdd=1
      nw=1
      nsz=neigd

      theta= 0.0
      phi= 0.0

      ALLOCATE ( ahelp(lmaxd*(lmaxd+2),natd,neigd,input%jspins) )
      ALLOCATE ( bhelp(lmaxd*(lmaxd+2),natd,neigd,input%jspins) )
      allocate ( hsomtx(neigd,neigd,2,2) )

      do ispin=1,input%jspins
       DO ie = 1, neigd
        DO na = 1, natd
         DO l = 1, lmaxd
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

      ALLOCATE(usdus%us(0:atoms%lmaxd,atoms%ntype,jspd), 
     +         usdus%dus(0:atoms%lmaxd,atoms%ntype,jspd),
     +         usdus%uds(0:atoms%lmaxd,atoms%ntype,jspd),
     +         usdus%duds(0:atoms%lmaxd,atoms%ntype,jspd),
     +         usdus%ddn(0:atoms%lmaxd,atoms%ntype,jspd),
     +         usdus%ulos(atoms%nlod,atoms%ntype,jspd),
     +         usdus%dulos(atoms%nlod,atoms%ntype,jspd),
     +         usdus%uulon(atoms%nlod,atoms%ntype,jspd),
     +         usdus%dulon(atoms%nlod,atoms%ntype,jspd))


      CALL spnorb(
     >         atoms,noco,input,mpi, enpara,
     >         vr,usdus,rsoc,.true.)
      
      do compo=1,3 
        rsoc%soangl= cmplx(0.0,0.0)

        
        CALL spnorb_angles(atoms,mpi,noco%theta,
     >        noco%phi,rsoc%soangl,compo)

        CALL hsoham(atoms,noco,input,nsz,neigd,chelp,rsoc,ahelp,
     >            bhelp,1,natd,mpi%n_rank,mpi%n_size,mpi%SUB_COMM,
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

      end subroutine wann_socmat_vec
      end module m_wann_socmat_vec