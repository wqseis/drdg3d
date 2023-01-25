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

module mod_plastic

  use mod_para,  only : RKIND, Np,                      &
                        cohesion, blkfric, Tvisc,       &
                        coef_bxx, coef_byy, coef_bxy,   &
                        fluidpres_profile_h1,           &
                        fluidpres_profile_h2,           &
                        fluidpres_profile_o1,           &
                        fluidpres_profile_o2
  use mod_types, only : meshvar
  use mod_eqns,  only : strain2stress, &
                        stress2strain

  implicit none

contains

!subroutine Return_Map(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs,depth,dt,&
!                      sxx,syy,szz,syz,sxz,sxy)
subroutine Return_Map(sxx,syy,szz,syz,sxz,sxy,rho,depth,dt)
  implicit none
  !real(kind=RKIND),intent(in) :: exx, eyy, ezz, exy, exz, eyz
  real(kind=RKIND),intent(inout) :: sxx, syy, szz, sxy, sxz, syz
  !real(kind=RKIND),intent(in) :: rho,cp,cs
  real(kind=RKIND),intent(in) :: rho, depth, dt

  real(kind=RKIND) :: sxx0, syy0, szz0, sxy0, sxz0, syz0, fluidpresh
  real(kind=RKIND) :: sm, sdxx, sdyy, sdzz, sdxy, sdxz, sdyz
  real(kind=RKIND) :: secinv, tau, taulim, decay, yldfac
  real(kind=RKIND) :: omeg, h1, h2, omeg1, omeg2
  real(kind=RKIND) :: byy, bxx, bxy

  !real(kind=RKIND) :: cohes, blkfric, angfric, tvisc!, dist
  real(kind=RKIND) :: angfric
  !integer :: i, j, k

  !real(kind=RKIND) :: lam,miu,chi

  ! For tpv27
  !cohes = 1.36e0
  !blkfric = 0.1934
  !tvisc = 0.03

  !!! For tpv29
  !cohes = 1.18e6
  !blkfric = 0.1680
  !tvisc = 0.05

  angfric = atan(blkfric)

  !do k = nk1, nk2
  !  do j = nj1, nj2
  !    do i = ni1, ni2
        !sxx = Txx(i,j,k); syy = Tyy(i,j,k); szz = Tzz(i,j,k)
        !sxy = Txy(i,j,k); sxz = Txz(i,j,k); syz = Tyz(i,j,k)
        !depth = -z(i,j,k)
        !depth = -z

  ! For tpv27
  !@fluidpresh = 1.0*9.8*depth ! rho*g*h
  !@if(depth<=15.0)then
  !@  omeg = 1.0
  !@elseif(depth<=20.0)then
  !@  omeg = (20.0 - depth)/5.0
  !@else
  !@  omeg = 0.0
  !@endif
  !@!szz0 = min(-2670.0*9.8*steph/3.0, -2670.0*9.8*depth)
  !@!szz0 = -2670.0*9.8*depth
  !@szz0 = -2.670*9.8*depth ! in MPa
  !@syy0 = omeg*(0.926793*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
  !@sxx0 = omeg*(1.073206*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
  !@sxy0 = omeg*(0.169029*(szz0+fluidpresh))
  !@sxz0 = 0.0
  !@syz0 = 0.0

  !h1 = 17
  !h2 = 22
  !omeg1 = 1.0
  !omeg2 = 0.0
  !byy = 1.025837
  !bxx = 0.974162
  !bxy = 0.158649
  !rho = 2.67
  h1 = fluidpres_profile_h1
  h2 = fluidpres_profile_h2
  omeg1 = fluidpres_profile_o1
  omeg2 = fluidpres_profile_o2
  byy = coef_byy
  bxx = coef_bxx
  bxy = coef_bxy

  ! For tpv30
  fluidpresh = 1.0*9.8*depth
  if(depth<=h1)then
    omeg = omeg1
  elseif(depth<=h2)then
    omeg = (h2-depth)/(h2-h1)
  else
    omeg = omeg2
  endif
  szz0 = min(-1.0e-3, -rho*9.8*depth) ! MPa
  syy0 = omeg*(byy*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
  sxx0 = omeg*(bxx*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
  sxy0 = omeg*(bxy*(szz0+fluidpresh))
  sxz0 = 0.0
  syz0 = 0.0

  !miu = rho*cs*cs
  !chi = rho*cp*cp
  !lam = chi-2d0*miu

  !sxx = sxx + chi*exx + lam*eyy + lam*ezz
  !syy = syy + lam*exx + chi*eyy + lam*ezz
  !szz = szz + lam*exx + lam*eyy + chi*ezz
  !sxy = sxy + miu*exy
  !sxz = sxz + miu*exz
  !syz = syz + miu*eyz

  sxx = sxx0 + sxx
  syy = syy0 + syy
  szz = szz0 + szz
  sxy = sxy0 + sxy
  sxz = sxz0 + sxz
  syz = syz0 + syz

  sm = (sxx+syy+szz)/3.0
  sdxx = sxx - sm
  sdyy = syy - sm
  sdzz = szz - sm
  sdxy = sxy
  sdxz = sxz
  sdyz = syz

  secinv = 0.5*(sdxx**2+sdyy**2+sdzz**2)+sdxy**2+sdxz**2+sdyz**2
  tau = sqrt(secinv)
  taulim = cohesion*cos(angfric) - (sm+fluidpresh)*sin(angfric)
  taulim = max(0.0, taulim)

  if(tau .gt. taulim) then
    decay = exp(-dt/tvisc)
    yldfac = decay + (1.0-decay)*taulim/tau
    sxx = sdxx*yldfac + sm
    syy = sdyy*yldfac + sm
    szz = sdzz*yldfac + sm
    sxy = sdxy*yldfac
    sxz = sdxz*yldfac
    syz = sdyz*yldfac
  endif

  sxx = sxx - sxx0
  syy = syy - syy0
  szz = szz - szz0
  sxy = sxy - sxy0
  sxz = sxz - sxz0
  syz = syz - syz0
  !Txx(i,j,k) = sxx; Tyy(i,j,k)= syy;  Tzz(i,j,k) = szz
  !Txy(i,j,k) = sxy; Txz(i,j,k)= sxz;  Tyz(i,j,k) = syz
  !    enddo
  !  enddo
  !enddo
endsubroutine

subroutine update_plastic(mesh,u)
  implicit none
  type(meshvar) :: mesh
  real(kind=RKIND),dimension(:,:) :: u ! (Np*Nelem,Nvar)
  real(kind=RKIND) :: dt,depth,rho,cp,cs
  real(kind=RKIND) :: exx,eyy,ezz,eyz,exz,exy
  real(kind=RKIND) :: sxx,syy,szz,syz,sxz,sxy
  integer :: ie,i

  dt = mesh%deltat

  do ie = 1,mesh%Nelem
    !do j = 1,Np
    !  iv(j) = j+(ie-1)*Np
    !end do
    cp = mesh%vp(ie)
    cs = mesh%vs(ie)
    rho = mesh%rho(ie)
    !rrho = 1d0/rho

    do i = 1,Np

      depth = -mesh%vz(i+(ie-1)*Np)
      exx = u(i+(ie-1)*Np,4)
      eyy = u(i+(ie-1)*Np,5)
      ezz = u(i+(ie-1)*Np,6)
      eyz = u(i+(ie-1)*Np,7)
      exz = u(i+(ie-1)*Np,8)
      exy = u(i+(ie-1)*Np,9)

      call strain2stress(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs, &
                         sxx,syy,szz,syz,sxz,sxy)

      call Return_Map(sxx,syy,szz,syz,sxz,sxy,rho,depth,dt)

      call stress2strain(sxx,syy,szz,syz,sxz,sxy,rho,cp,cs, &
                         exx,eyy,ezz,eyz,exz,exy)

      u(i+(ie-1)*Np,4) = exx
      u(i+(ie-1)*Np,5) = eyy
      u(i+(ie-1)*Np,6) = ezz
      u(i+(ie-1)*Np,7) = eyz
      u(i+(ie-1)*Np,8) = exz
      u(i+(ie-1)*Np,9) = exy
    end do
  end do ! element

end subroutine

end module
