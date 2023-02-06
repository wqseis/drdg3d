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

module mod_interp

  use mod_para,   only : RKIND
  use mod_math,   only : cross

  implicit none

contains

function tri_interp(v1,v2,v3,c1,c2,c3,p) result (cp)

  real(kind=rkind),dimension(3) :: v1,v2,v3,p
  real(kind=rkind) :: c1,c2,c3,cp
  cp = tri_interp_barycentric(v1,v2,v3,c1,c2,c3,p)

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

function tri_interp_barycentric(v1,v2,v3,c1,c2,c3,p) result (cp)
  implicit none

  real(kind=rkind),dimension(3) :: v1,v2,v3
  real(kind=rkind) :: c1,c2,c3,cp
  real(kind=rkind),dimension(3) :: p
  real(kind=rkind),dimension(3) :: n
  real(kind=rkind) :: w1,w2,w3
  real(kind=rkind) :: Xv1,Xv2,Xv3
  real(kind=rkind) :: Yv1,Yv2,Yv3
  real(kind=rkind) :: Px,Py,norm

  n = cross(v2-v1,v3-v1)
  if (abs(n(1))>abs(n(2)) .and. abs(n(1))>abs(n(3))) then
    Xv1=v1(2)
    Xv2=v2(2)
    Xv3=v3(2)
    Yv1=v1(3)
    Yv2=v2(3)
    Yv3=v3(3)
    Px=p(2)
    Py=p(3)
  else if( abs(n(2))>abs(n(1)) .and. abs(n(2))>abs(n(3))) then
    Xv1=v1(1)
    Xv2=v2(1)
    Xv3=v3(1)
    Yv1=v1(3)
    Yv2=v2(3)
    Yv3=v3(3)
    Px=p(1)
    Py=p(3)
  else
    Xv1=v1(1)
    Xv2=v2(1)
    Xv3=v3(1)
    Yv1=v1(2)
    Yv2=v2(2)
    Yv3=v3(2)
    Px=p(1)
    Py=p(2)
  end if

  ! barycentric
  norm=(Yv2-Yv3)*(Xv1-Xv3)+(Xv3-Xv2)*(Yv1-Yv3)
  w1=  (Yv2-Yv3)*(Px -Xv3)+(Xv3-Xv2)*(Py -Yv3)
  w2=  (Yv3-Yv1)*(Px -Xv3)+(Xv1-Xv3)*(Py -Yv3)
  w1=w1/norm
  w2=w2/norm
  w3=1d0-w1-w2

  cp=(w1*c1+w2*c2+w3*c3)/(w1+w2+w3)

end function

function interp_dist_n(x,y,z,v,n,x1,y1,z1) result (v1)
  implicit none

  integer :: n,i
  real(kind=rkind),dimension(n) :: x,y,z,v
  real(kind=rkind) :: x1,y1,z1,v1
  real(kind=rkind) :: r1!,r2,r3
  real(kind=rkind) :: w1!,w2,w3,
  real(kind=rkind) :: wsum,vsum

  vsum=0.
  wsum=0.
  do i = 1,n
    r1=sqrt((x(i)-x1)**2+(y(i)-y1)**2+(z(i)-z1)**2)
    w1=1d0/(r1+1d-30)
    vsum=vsum+w1*v(i)
    wsum=wsum+w1
  end do

  v1=vsum/wsum

end function


end module
