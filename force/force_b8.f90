!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_forceb8
   !-----------------------------------------------------------------------------
   ! Implements the surface contribution to the force. Madsen Eq.(B8)
   !
   ! FZJ 15/3-01 GMadsen
   !-----------------------------------------------------------------------------
CONTAINS
   SUBROUTINE force_b8(atoms,ecwk,stars, sym,cell, jspin, force,f_b8)

      USE m_constants, ONLY : tpi_const
      USE m_sphbes
      USE m_stern
      USE m_types

      IMPLICIT NONE

      TYPE(t_sym),   INTENT(IN) :: sym
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_atoms), INTENT(IN) :: atoms

      INTEGER, INTENT(IN)    :: jspin
      COMPLEX, INTENT(IN)    :: ecwk(stars%ng3)
      COMPLEX, INTENT(INOUT) :: f_b8(3,atoms%ntype)
      REAL,    INTENT(INOUT) :: force(:,:,:)

      INTEGER g(3),nst,stg(3,sym%nop),ia,istr,i,j,jj,iType
      REAL    fj(0:atoms%lmaxd),rotkzz(3),rstg(3,sym%nop)
      REAL    frmt,gl,pha,s
      COMPLEX taup(sym%nop),factor,fact,fstar(3),fsur2(3)
   
      DO iType=1,atoms%ntype
         ia = atoms%firstAtom(iType)
         frmt = 2.0*tpi_const*atoms%rmt(iType)**2
         fsur2(1:3) = cmplx(0.0,0.0)         

         ! Skip G=(0,0,0) [no contribution to ekin]
         DO istr=2,stars%ng3_fft 
            g(:)     = stars%kv3(:,istr)
            fstar(:) = cmplx(0.0,0.0)
            CALL stern(sym,cell,g, nst,stg,taup,gl,rstg)

            CALL sphbes(atoms%lmax(iType),atoms%rmt(iType)*gl,fj)
            fact = ecwk(istr) * fj(1) / gl

            DO jj=1,nst
               pha=(atoms%taual(1,ia)*stg(1,jj)+atoms%taual(2,ia)*stg(2,jj)&
                   +atoms%taual(3,ia)*stg(3,jj))*tpi_const

               ! Swapped sin and cos because there's an i in the equation
               factor = fact * cmplx(-sin(pha),cos(pha)) * taup(jj)
               DO i=1,3
                  fstar(i) = fstar(i) + factor*rstg(i,jj)
               END DO
            END DO
            DO i=1,3
               fsur2(i)=fsur2(i)+fstar(i)*frmt
            END DO
         END DO
         DO i=1,3
            f_b8(i,iType) = f_b8(i,iType) + fsur2(i)
            force(i,iType,jspin) = force(i,iType,jspin) + real(fsur2(i))
         END DO
      END DO

   END SUBROUTINE force_b8
END MODULE m_forceb8
