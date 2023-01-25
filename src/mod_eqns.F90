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

module mod_eqns
  use mod_para, only : RKIND
  implicit none

contains

!function norm3(a)
!  real(kind=RKIND) :: a(3)
!  real(kind=RKIND) :: norm3
!  !norm3=dsqrt(a(1)**2+a(2)**2+a(3)**2)
!  norm3=sqrt(a(1)**2+a(2)**2+a(3)**2)
!end

subroutine Flux1(U,rho,cp,cs,F)
  implicit none
  real(kind=rkind),intent(in) :: U(:),rho,cp,cs
  real(kind=rkind),intent(out) :: F(:)
  real(kind=rkind) :: rrho,miu,lam,chi
  rrho = 1d0/rho
  miu = rho*cs**2
  chi = rho*cp**2
  lam = chi-2d0*miu

  ! U=(rho*Vx,rho*Vy,rho*Vz,Exx,Eyy,Ezz,Eyz,Exz,Exy)
  ! F=(Sxx,Sxy,Sxz,Vx,0,0,0,Vz,Vy)
  F(1) = U(4) * chi + U(5) * lam + U(6) * lam
  F(2) = U(9) * miu
  F(3) = U(8) * miu
  F(4) = U(1) * rrho
  F(5) = 0
  F(6) = 0
  F(7) = 0
  F(8) = U(3) * rrho
  F(9) = U(2) * rrho
end subroutine

subroutine Flux2(U,rho,cp,cs,F)
  implicit none
  real(kind=rkind),intent(in) :: U(:),rho,cp,cs
  real(kind=rkind),intent(out) :: F(:)
  real(kind=rkind) :: rrho,miu,lam,chi
  rrho = 1d0/rho
  miu = rho*cs**2
  chi = rho*cp**2
  lam = chi-2d0*miu

  ! U=(rho*Vx,rho*Vy,rho*Vz,Exx,Eyy,Ezz,Eyz,Exz,Exy)
  ! F=(Sxy,Syy,Syz,0,Vy,0,Vz,0,Vx)
  F(1) = U(9) * miu
  F(2) = U(4) * lam + U(5) * chi + U(6) * lam
  F(3) = U(7) * miu
  F(4) = 0
  F(5) = U(2) * rrho
  F(6) = 0
  F(7) = U(3) * rrho
  F(8) = 0
  F(9) = U(1) * rrho
end subroutine

subroutine Flux3(U,rho,cp,cs,F)
  implicit none
  real(kind=rkind),intent(in) :: U(:),rho,cp,cs
  real(kind=rkind),intent(out) :: F(:)
  real(kind=rkind) :: rrho,miu,lam,chi
  rrho = 1d0/rho
  miu = rho*cs**2
  chi = rho*cp**2
  lam = chi-2d0*miu

  ! U=(rho*Vx,rho*Vy,rho*Vz,Exx,Eyy,Ezz,Eyz,Exz,Exy)
  ! F=(Sxz,Syz,Szz,0,0,Vz,Vy,Vx,0)
  F(1) = U(8) * miu
  F(2) = U(7) * miu
  F(3) = U(4) * lam + U(5) * lam + U(6) * chi
  F(4) = 0
  F(5) = 0
  F(6) = U(3) * rrho
  F(7) = U(2) * rrho
  F(8) = U(1) * rrho
  F(9) = 0
end subroutine

subroutine strain2stress(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs,&
                         sxx,syy,szz,syz,sxz,sxy)
  implicit none
  real(kind=RKIND),intent(in) :: exx,eyy,ezz,eyz,exz,exy,rho,cp,cs
  real(kind=RKIND),intent(out) :: sxx,syy,szz,syz,sxz,sxy
  !real(kind=rkind) :: rrho,miu,lam,chi
  real(kind=rkind) :: miu,lam,chi
  !rrho = 1d0/rho
  miu = rho*cs**2
  chi = rho*cp**2
  lam = chi-2d0*miu

  sxx = exx * chi + eyy * lam + ezz * lam
  syy = exx * lam + eyy * chi + ezz * lam
  szz = exx * lam + eyy * lam + ezz * chi
  syz = eyz * miu
  sxz = exz * miu
  sxy = exy * miu
end subroutine

subroutine stress2strain(sxx,syy,szz,syz,sxz,sxy,rho,cp,cs,&
                         exx,eyy,ezz,eyz,exz,exy)
  implicit none
  real(kind=RKIND),intent(in) :: sxx,syy,szz,syz,sxz,sxy,rho,cp,cs
  real(kind=RKIND),intent(out) :: exx,eyy,ezz,eyz,exz,exy
  !real(kind=rkind) :: rrho,miu,lam,chi
  real(kind=rkind) :: miu,lam,chi,a,b
  !rrho = 1d0/rho
  miu = rho*cs**2
  chi = rho*cp**2
  lam = chi-2d0*miu
  a = (chi+lam)/(chi**2+lam*chi-2.0*lam**2)
  b = (-lam   )/(chi**2+lam*chi-2.0*lam**2)

  exx = sxx * a + syy * b + szz * b
  eyy = sxx * b + syy * a + szz * b
  ezz = sxx * b + syy * b + szz * a
  eyz = syz / miu
  exz = sxz / miu
  exy = sxy / miu
end subroutine

subroutine generate_fluctuations(z,T,T_hat,v,v_hat,F)
  implicit none
  real(kind=rkind),intent(in) :: z,T,T_hat,v,v_hat
  real(kind=rkind),intent(out) :: F
  F = 0.5d0*(z*(v_hat-v) + (T_hat-T))
end subroutine

subroutine extract_traction_velocity(u,n,rho,cp,cs,vx,vy,vz,Tx,Ty,Tz)
  implicit none
  real(kind=rkind),intent(in) :: u(9),n(3),rho,cp,cs
  real(kind=rkind),intent(out) :: vx,vy,vz,Tx,Ty,Tz
  real(kind=rkind) :: rrho
  real(kind=rkind) :: miu,lam,chi
  real(kind=rkind) :: exx,eyy,ezz,eyz,exz,exy
  real(kind=rkind) :: sxx,syy,szz,syz,sxz,sxy

  rrho = 1.0/rho
  miu = rho*cs**2
  chi = rho*cp**2
  lam = chi-2d0*miu

  vx  = u(1)*rrho
  vy  = u(2)*rrho
  vz  = u(3)*rrho
  exx = u(4)
  eyy = u(5)
  ezz = u(6)
  eyz = u(7)
  exz = u(8)
  exy = u(9)

  sxx = exx * chi + eyy * lam + ezz * lam
  syy = exx * lam + eyy * chi + ezz * lam
  szz = exx * lam + eyy * lam + ezz * chi
  syz = eyz * miu
  sxz = exz * miu
  sxy = exy * miu

  Tx = sxx*n(1) + sxy*n(2) + sxz*n(3)
  Ty = sxy*n(1) + syy*n(2) + syz*n(3)
  Tz = sxz*n(1) + syz*n(2) + szz*n(3)
end subroutine

subroutine riemannSolver_continuous(v_p,v_m,sigma_p,sigma_m,z_p,z_m, &
    v_hat_p,v_hat_m,sigma_hat_p,sigma_hat_m)
  implicit none
  real(kind=rkind),intent(in) :: v_p,v_m,sigma_p,sigma_m,z_p,z_m
  real(kind=rkind),intent(out) :: v_hat_p,v_hat_m,sigma_hat_p,sigma_hat_m
  real(kind=rkind) :: p
  real(kind=rkind) :: q
  real(kind=rkind) :: phi
  real(kind=rkind) :: eta

  p = z_m*v_p + sigma_p
  q = z_p*v_m - sigma_m

  eta = (z_p*z_m)/(z_p+z_m)

  phi = eta*(p/z_p - q/z_m)

  sigma_hat_p = phi
  sigma_hat_m = phi

  v_hat_p = (q+phi)/z_m
  v_hat_m = (p-phi)/z_p
end subroutine

end module
