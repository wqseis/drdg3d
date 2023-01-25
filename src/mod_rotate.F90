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

module mod_rotate

  use mod_para,only : RKIND
  implicit none

contains

subroutine rotation_matrix_velocity(n,s,t,Rot)
  implicit none
  real(kind=RKIND),intent(in) :: n(3),s(3),t(3)
  real(kind=RKIND),intent(out) :: Rot(3,3)
  Rot(1,:) = n
  Rot(2,:) = s
  Rot(3,:) = t
end subroutine

subroutine rotation_matrix_velocity_inv(n,s,t,Rot)
  implicit none
  real(kind=RKIND),intent(in) :: n(3),s(3),t(3)
  real(kind=RKIND),intent(out) :: Rot(3,3)
  Rot(:,1) = n
  Rot(:,2) = s
  Rot(:,3) = t
end subroutine

subroutine rotate_xyz2nml(n,m,l,Tx,Ty,Tz,Tn,Tm,Tl)
  implicit none
  real(kind=RKIND),intent(in) :: n(3),m(3),l(3),Tx,Ty,Tz
  real(kind=RKIND),intent(out) :: Tn,Tm,Tl
  Tn = Tx*n(1)+Ty*n(2)+Tz*n(3)
  Tm = Tx*m(1)+Ty*m(2)+Tz*m(3)
  Tl = Tx*l(1)+Ty*l(2)+Tz*l(3)
end subroutine

subroutine rotate_nml2xyz(n,m,l,Tn,Tm,Tl,Tx,Ty,Tz)
  implicit none
  real(kind=RKIND),intent(in) :: n(3),m(3),l(3), Tn,Tm,Tl
  real(kind=RKIND),intent(out) :: Tx,Ty,Tz
  Tx = Tn*n(1)+Tm*m(1)+Tl*l(1)
  Ty = Tn*n(2)+Tm*m(2)+Tl*l(2)
  Tz = Tn*n(3)+Tm*m(3)+Tl*l(3)
end subroutine

subroutine rotate_u(Tv,Ts,u)
  implicit none
  real(kind=RKIND),intent(inout) :: u(9)
  real(kind=RKIND),intent(in) :: Ts(6,6),Tv(3,3)
  real(kind=RKIND) :: u1(9)

  !u(1:3) = matmul(Tv,u(1:3))
  u1(1) = Tv(1,1)*u(1)+Tv(1,2)*u(2)+Tv(1,3)*u(3)
  u1(2) = Tv(2,1)*u(1)+Tv(2,2)*u(2)+Tv(2,3)*u(3)
  u1(3) = Tv(3,1)*u(1)+Tv(3,2)*u(2)+Tv(3,3)*u(3)
  !u(4:9) = matmul(Ts,u(4:9))
  u1(4) = Ts(1,1)*u(4)+Ts(1,2)*u(5)+Ts(1,3)*u(6)+Ts(1,4)*u(7)+Ts(1,5)*u(8)+Ts(1,6)*u(9)
  u1(5) = Ts(2,1)*u(4)+Ts(2,2)*u(5)+Ts(2,3)*u(6)+Ts(2,4)*u(7)+Ts(2,5)*u(8)+Ts(2,6)*u(9)
  u1(6) = Ts(3,1)*u(4)+Ts(3,2)*u(5)+Ts(3,3)*u(6)+Ts(3,4)*u(7)+Ts(3,5)*u(8)+Ts(3,6)*u(9)
  u1(7) = Ts(4,1)*u(4)+Ts(4,2)*u(5)+Ts(4,3)*u(6)+Ts(4,4)*u(7)+Ts(4,5)*u(8)+Ts(4,6)*u(9)
  u1(8) = Ts(5,1)*u(4)+Ts(5,2)*u(5)+Ts(5,3)*u(6)+Ts(5,4)*u(7)+Ts(5,5)*u(8)+Ts(5,6)*u(9)
  u1(9) = Ts(6,1)*u(4)+Ts(6,2)*u(5)+Ts(6,3)*u(6)+Ts(6,4)*u(7)+Ts(6,5)*u(8)+Ts(6,6)*u(9)

  u = u1
end subroutine

subroutine rotation_matrix_strain(n,s,t,Rot)
  implicit none
  real(kind=RKIND),intent(in) :: n(3),s(3),t(3)
  real(kind=RKIND),intent(out) :: Rot(6,6)
  real(kind=RKIND) :: nx,ny,nz,sx,sy,sz,tx,ty,tz
  nx = n(1); ny = n(2); nz = n(3)
  sx = s(1); sy = s(2); sz = s(3)
  tx = t(1); ty = t(2); tz = t(3)
  Rot(1,:) = (/  nx*nx,   ny*ny,   nz*nz,         ny*nz,         nx*nz,         nx*ny/)
  Rot(2,:) = (/  sx*sx,   sy*sy,   sz*sz,         sy*sz,         sx*sz,         sx*sy/)
  Rot(3,:) = (/  tx*tx,   ty*ty,   tz*tz,         ty*tz,         tx*tz,         tx*ty/)
  Rot(4,:) = (/2*sx*tx, 2*sy*ty, 2*sz*tz, sy*tz + sz*ty, sx*tz + sz*tx, sx*ty + sy*tx/)
  Rot(5,:) = (/2*nx*tx, 2*ny*ty, 2*nz*tz, ny*tz + nz*ty, nx*tz + nz*tx, nx*ty + ny*tx/)
  Rot(6,:) = (/2*nx*sx, 2*ny*sy, 2*nz*sz, ny*sz + nz*sy, nx*sz + nz*sx, nx*sy + ny*sx/)
end subroutine

subroutine rotation_matrix_strain_inv(n,s,t,Rot)
  implicit none
  real(kind=RKIND),intent(in) :: n(3),s(3),t(3)
  real(kind=RKIND),intent(out) :: Rot(6,6)
  real(kind=RKIND) :: nx,ny,nz,sx,sy,sz,tx,ty,tz
  nx = n(1); ny = n(2); nz = n(3)
  sx = s(1); sy = s(2); sz = s(3)
  tx = t(1); ty = t(2); tz = t(3)
  Rot(1,:) = (/  nx*nx,   sx*sx,   tx*tx,         sx*tx,         nx*tx,         nx*sx/)
  Rot(2,:) = (/  ny*ny,   sy*sy,   ty*ty,         sy*ty,         ny*ty,         ny*sy/)
  Rot(3,:) = (/  nz*nz,   sz*sz,   tz*tz,         sz*tz,         nz*tz,         nz*sz/)
  Rot(4,:) = (/2*ny*nz, 2*sy*sz, 2*ty*tz, sy*tz + sz*ty, ny*tz + nz*ty, ny*sz + nz*sy/)
  Rot(5,:) = (/2*nx*nz, 2*sx*sz, 2*tx*tz, sx*tz + sz*tx, nx*tz + nz*tx, nx*sz + nz*sx/)
  Rot(6,:) = (/2*nx*ny, 2*sx*sy, 2*tx*ty, sx*ty + sy*tx, nx*ty + ny*tx, nx*sy + ny*sx/)
end subroutine

subroutine inverse_rotation_matrix_strain(R,invR)
  implicit none
  real(kind=RKIND),intent(in) :: R(6,6)
  real(kind=RKIND),intent(out) :: invR(6,6)
  invR = transpose(R);
  !invR(1:3,4:6) = 2.0*invR(1:3,4:6);
  !invR(4:6,1:3) = 0.5*invR(4:6,1:3);
  invR(1:3,4:6) = 0.5*invR(1:3,4:6);
  invR(4:6,1:3) = 2.0*invR(4:6,1:3);
end subroutine

subroutine inverse_rotation_matrix_velocity(R,invR)
  implicit none
  real(kind=RKIND),intent(in) :: R(3,3)
  real(kind=RKIND),intent(out) :: invR(3,3)
  invR = transpose(R);
end subroutine

end module
