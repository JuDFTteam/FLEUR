!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_eulerrot
c************************************
c     Perform Euler rotations.
c     Y. Mokrousov 
c
c     tidied up version
c     Frank Freimuth
c************************************
      contains
      subroutine eulerrot(
     >               nwfs,alpha,beta,gamma,
     <               amx)
c************************************************************************
c..Perform "nwfs" Euler rotations.
c..Input: Set of "nwfs" Euler angles: alpha(i),beta(i),gamma(i); i=1,nwfs
c..Output: "nwfs" rotation matrices: amx(:,:,i); i=1,nwfs
c************************************************************************

      implicit none

      integer, intent (in)  :: nwfs
      real,    intent (in)  :: alpha(nwfs),beta(nwfs),gamma(nwfs) 
      real,    intent (out) :: amx(3,3,nwfs)

      integer :: nwf

      do nwf = 1,nwfs
         call eulerrot1(
     >                  alpha(nwf),beta(nwf),gamma(nwf),
     <                  amx(1,1,nwf) )
      enddo 

      end subroutine eulerrot

      subroutine eulerrot1(
     >               alpha,beta,gamma,
     <               amx)
c***********************************************************************
c..Perform one Euler rotation.
c..Input: Euler angles: alpha, beta, gamma
c..Output: Rotation matrix: amx(:,:)
c 
c..Given the Euler angles, the following procedure is 
c..performed:
c..1. Rotation around the z-axis by alpha (matrix D)
c..2. Rotation around the x-axis by beta  (matrix C)
c..3. Rotation around the z-axis by gamma again. (matrix B)
c..The overall rotation is given by the matrix A = B*C*D
c***********************************************************************

      implicit none

      real,    intent (in)  :: alpha,beta,gamma 
      real,    intent (out) :: amx(3,3)

      real    :: bmx(3,3),cmx(3,3),dmx(3,3),hmx(3,3)

      dmx(1,1) = cos(alpha) ; dmx(1,2) = sin(alpha) ; dmx(1,3) = 0. 
      dmx(2,1) =-sin(alpha) ; dmx(2,2) = cos(alpha) ; dmx(2,3) = 0. 
      dmx(3,1) = 0.         ; dmx(3,2) = 0.         ; dmx(3,3) = 1. 

      cmx(1,1) = 1.  ; cmx(1,2) = 0.         ; cmx(1,3) = 0. 
      cmx(2,1) = 0.  ; cmx(2,2) = cos(beta)  ; cmx(2,3) = sin(beta)
      cmx(3,1) = 0.  ; cmx(3,2) =-sin(beta)  ; cmx(3,3) = cos(beta)
  
      bmx(1,1) = cos(gamma) ; bmx(1,2) = sin(gamma) ; bmx(1,3) = 0. 
      bmx(2,1) =-sin(gamma) ; bmx(2,2) = cos(gamma) ; bmx(2,3) = 0. 
      bmx(3,1) = 0.         ; bmx(3,2) = 0.         ; bmx(3,3) = 1. 

      hmx = matmul(cmx,dmx)
      amx = matmul(bmx,hmx)

      end subroutine eulerrot1


      end module m_eulerrot  
