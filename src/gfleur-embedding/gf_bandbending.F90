!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_gf_bandbending
!Created by wortmann on May 9, 2012
      implicit none
      contains

      subroutine gf_bandbending(lapw,jspins,nkpts)
      use m_gf_types
      use m_gf_energies
      use m_gf_io2dmat
      use m_gf_iotmat
      USE m_gf_propaemb
      use m_gf_embedding
      USE m_gf_math
      USE m_gf_phasematrix
      implicit none
      integer,intent(in)      :: jspins,nkpts
      type(t_lapw),intent(in) :: lapw

      logical :: l_exist
      integer :: no_layers,n,nk,en,en_mod,jspin
      real,allocatable :: bend(:),energies(:)
      complex,allocatable :: embpot_in(:,:),embpot_out(:,:)
      complex,allocatable :: tmat(:,:)



      inquire(file="gf_bandbending",exist=l_exist)
      if (.not.l_exist) return


      allocate(energies(gf_noen()))
      allocate(embpot_in(lapw%nv2_tot,lapw%nv2_tot))
      allocate(embpot_out(lapw%nv2_tot,lapw%nv2_tot))
      allocate(tmat(2*lapw%nv2_tot,2*lapw%nv2_tot))
      energies=real(gf_allz(1))


      write(*,*) "Apply bandbending to embpot"
      open(999,file="gf_bandbending")
      read(999,*) no_layers
      allocate(bend(no_layers))
      DO n=1,no_layers
          read(999,*) bend(n)
      ENDDO
      close(999)


      DO jspin=1,jspins
        DO nk=1,nkpts
            DO en=1,gf_noen()
                en_mod=get_energy_index(energies(en)-bend(1),energies)
                IF (en_mod==-1) THEN
                     write(*,*) "Energy out of range:",energies(en)
                     CYCLE
                ENDIF
                CALL gf_GETEMB2(embpot_in,1,1,en_mod,nk,jspin,lapw)
                DO n=1,no_layers
                    en_mod=get_energy_index(energies(en)-bend(n),energies)
                    CALL gf_read_tmat2(                                               &
     &                   1,en_mod,nk,jspin,lapw,                                          &
     &                   tmat)
                    CALL gf_propaemb(.FALSE.,lapw%nv2_tot,embpot_in,   &
     &                      tmat,embpot_out)
                    embpot_in=matmul(getPhaseMatrix(),embpot_in)
                    embpot_out=matmul(embpot_in,mat_inverse(getPhaseMatrix()))
                ENDDO
                CALL gf_write2Dmat(IO2D_EMB,2,1,en,nk,jspin,&
     &               lapw,embpot_in)
           ENDDO
        ENDDO
      ENDDO

      CONTAINS

      FUNCTION get_energy_index(z,energies)
      implicit none
      real,intent(in):: z,energies(:)
      integer        :: get_energy_index

      integer :: n(size(energies))
      get_energy_index=-1

      if (all(energies<z).or.all(energies>z)) return

      n=minloc(abs(z-energies))

      get_energy_index=n(1)

      end function

    end subroutine



end module m_gf_bandbending
