!-----------------------------------------------------------------------
!Copyright (C) 2021-2023 Wenqiang ZHANG (wqseis@gmail.com)
!
!This file is part of DRDG3D.
!
!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------

module mod_math

  use mod_para, only : RKIND

  implicit none

contains

!function norm3(a)
!  real(kind=RKIND) :: a(3)
!  real(kind=RKIND) :: norm3
!  !norm3=dsqrt(a(1)**2+a(2)**2+a(3)**2)
!  norm3=sqrt(a(1)**2+a(2)**2+a(3)**2)
!end

function cross(a, b)
  implicit none
  real(kind=RKIND), dimension(3) :: cross
  real(kind=RKIND), dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
end function

function tri_interp_dist(v1,v2,v3,c1,c2,c3,p) result (cp)
  implicit none

  real(kind=rkind),dimension(3) :: v1,v2,v3
  real(kind=rkind) :: c1,c2,c3,cp
  real(kind=rkind),dimension(3) :: p
  real(kind=rkind) :: r1,r2,r3
  real(kind=rkind) :: w1,w2,w3

  r1=sqrt((p(1)-v1(1))**2+(p(2)-v1(2))**2+(p(3)-v1(3))**2)
  r2=sqrt((p(1)-v2(1))**2+(p(2)-v2(2))**2+(p(3)-v2(3))**2)
  r3=sqrt((p(1)-v3(1))**2+(p(2)-v3(2))**2+(p(3)-v3(3))**2)

  w1=1d0/(r1+1d-300)
  w2=1d0/(r2+1d-300)
  w3=1d0/(r3+1d-300)

  cp=(w1*c1+w2*c2+w3*c3)/(w1+w2+w3)

end function

subroutine mxm(a,b,c,m,n,l)
  implicit none
  integer :: i,j,k,m,n,l
  real(kind=rkind),dimension(:,:) :: a,b,c
  c=0.
  do i=1,m
    do j=1,l
      do k=1,n
        c(i,j) = c(i,j) + a(i,k)*b(k,j)
      end do
    end do
  end do
end subroutine

end module
