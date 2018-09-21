!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_socmat
      contains
      subroutine wann_socmat(
     >               mpi,enpara,input,noco,atoms,
     >               lmaxd,ntypd,nlod,natd,neigd,
     >               llod,jmtd,jspd,nlhd,neq,
     >               ntype,theta_in,phi_in,jspins,irank,
     >               jri,lmax,dx,rmsh,el,ello,nlo,llo,
     >               l_dulo,ulo_der,vr,
     >               acof,bcof,chelp,
     <               hsomtx)
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
      integer, intent(in) :: ntypd
      integer, intent(in) :: nlod
      integer, intent(in) :: natd
      integer, intent(in) :: neigd

      integer, intent(in) :: llod
      integer, intent(in) :: jmtd
      integer, intent(in) :: jspd
      integer, intent(in) :: nlhd
      integer, intent(in) :: neq(:) !ntype

      integer, intent(in) :: ntype
      real,    intent(in) :: theta_in
      real,    intent(in) :: phi_in
      integer, intent(in) :: jspins
      integer, intent(in) :: irank

      integer, intent(in) :: jri(:)  ! ntypd
      integer, intent(in) :: lmax(:) ! ntypd
      real,    intent(in) :: dx(:)   ! ntypd
      real,    intent(in) :: rmsh(:,:)  ! jmtd,ntypd
      real,    intent(in) :: el(0:,:,:) ! 0:lmaxd,ntypd,max(2,jspd)
      real,    intent(in) :: ello(:,:,:)! nlod,ntypd,max(2,jspd)
      integer, intent(in) :: nlo(:)     ! ntypd
      integer, intent(in) :: llo(:,:)   ! nlod,ntypd

      logical, intent(in) :: l_dulo(:,:)  ! l_dulo(nlod,ntypd)
      integer, intent(in) :: ulo_der(:,:) ! ulo_der(nlod,ntypd)
      real,    intent(in) :: vr(:,0:,:,:)    ! jmtd,0:nlhd,ntypd,jspd   
  
      complex, intent(in) :: acof(:,0:,:,:) !acof(noccbd,0:lmd,natd,jspd)
      complex, intent(in) :: bcof(:,0:,:,:) !bcof(noccbd,0:lmd,natd,jspd)
      complex, intent(in) :: chelp(-llod:,:,:,:,:) !chelp(-llod:llod,neigd,nlod,natd,jspd)

      complex, allocatable, intent(out):: hsomtx(:,:,:,:) !(neigd,neigd,2,2)

      integer :: n,l,nwdd,nw,ispin,ie,na,ll1,m,lm,i,nsz(2)
      real    :: s(3),r2
      logical :: l_all
      CHARACTER*3 chntype
      real    :: theta,phi,pi

      real,allocatable :: ddn(:,:,:) ! 0:lmaxd,ntypd,jspd
      real,allocatable :: us(:,:,:)  ! 0:lmaxd,ntypd,jspd
      real,allocatable :: dus(:,:,:) ! 0:lmaxd,ntypd,jspd
      real,allocatable :: uds(:,:,:) ! 0:lmaxd,ntypd,jspd
      real,allocatable :: duds(:,:,:)! 0:lmaxd,ntypd,jspd

      real,allocatable :: ulos(:,:,:)  ! nlod,ntypd,jspd
      real,allocatable :: dulos(:,:,:) ! nlod,ntypd,jspd
      real,allocatable :: uulon(:,:,:) ! nlod,ntypd,jspd
      real,allocatable :: dulon(:,:,:) ! nlod,ntypd,jspd  

      COMPLEX,ALLOCATABLE :: ahelp(:,:,:,:)
      complex,allocatable :: bhelp(:,:,:,:)

      TYPE(t_rsoc)::rsoc
      TYPE(t_usdus):: usdus

      nwdd=1
      nw=1
      nsz=neigd

      IF (.true.) THEN
        theta= theta_In
        phi= phi_In
      ELSE
        pi= 4.*ATAN(1.)
        theta= -theta_In
        phi=   phi_In+pi
        ! now the definition of rotation matrices
        ! is equivalent to the def in the noco-routines
      ENDIF

      ALLOCATE ( ahelp(lmaxd*(lmaxd+2),natd,neigd,input%jspins) )
      ALLOCATE ( bhelp(lmaxd*(lmaxd+2),natd,neigd,input%jspins) )

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

      ALLOCATE( us(0:lmaxd,ntypd,jspd), dus(0:lmaxd,ntypd,jspd),
     +          uds(0:lmaxd,ntypd,jspd),duds(0:lmaxd,ntypd,jspd),
     +          ddn(0:lmaxd,ntypd,jspd),
     +          ulos(nlod,ntypd,jspd),dulos(nlod,ntypd,jspd),
     +          uulon(nlod,ntypd,jspd),dulon(nlod,ntypd,jspd))

      CALL spnorb(
     >         atoms,noco,input,mpi, enpara,
     >         vr,usdus,rsoc,.true.)
      
      rsoc%soangl= conjg(rsoc%soangl)

      CALL hsoham(atoms,noco,input,nsz,neigd,chelp,rsoc,ahelp,
     >            bhelp,1,natd,mpi%n_rank,mpi%n_size,mpi%SUB_COMM,
     <            hsomtx)
 
      end subroutine wann_socmat
      end module m_wann_socmat
