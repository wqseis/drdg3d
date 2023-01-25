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

module mod_rk

use mod_para,  only : RKIND,               &
                      rk2a, rk2b,          &
                      rk3a, rk3b,          &
                      rk4a, rk4b,          &
                      timeIntegrationMethod
use mod_types, only : meshvar

implicit none

contains


subroutine rk_init(mesh)
  implicit none
  type(meshvar) :: mesh
  ! RK method
  !mesh%nrk = 4
  if(timeIntegrationMethod == 0) then
    mesh%nrk = 4
  elseif(timeIntegrationMethod == 1) then
    mesh%nrk = 5
  elseif(timeIntegrationMethod == 2) then
    mesh%nrk = 3
  elseif(timeIntegrationMethod == 3) then
    mesh%nrk = 2
  end if
end subroutine

!    mu = u
!    mesh%mslip = mesh%slip
!    tu = 0.0
!    mesh%tslip = 0.0
!    ! rate state
!    mesh%mstate = mesh%state
!    mesh%tstate = 0.0

subroutine rk_update(mesh,u,mu,tu,hu,dt,irk)
  implicit none

  type(meshvar) :: mesh
  real(kind=rkind),dimension(:,:) :: u,mu,tu,hu

  integer :: irk
  real(kind=rkind) :: dt

  if (timeIntegrationMethod == 0) then
    if(irk==1)then
      u  = mu+0.5d0*dt*hu
      tu = mu+1.0d0/6.0d0*dt*hu

      !s  = ms+0.5*dt*hs
      !ts = ms+1.0/6.0*dt*hs
      mesh%slip  = mesh%mslip+0.5d0*dt*mesh%sliprate
      mesh%tslip = mesh%mslip+1.0d0/6.0d0*dt*mesh%sliprate
      ! rate state
      mesh%state  = mesh%mstate+0.5d0*dt*mesh%hstate
      mesh%tstate = mesh%mstate+1.0d0/6.0d0*dt*mesh%hstate
    elseif(irk==2)then
      u  = mu+0.5d0*dt*hu
      tu = tu+1.0d0/3.0d0*dt*hu

      mesh%slip  = mesh%mslip+0.5d0*dt*mesh%sliprate
      mesh%tslip = mesh%tslip+1.0d0/3.0d0*dt*mesh%sliprate
      ! rate state
      mesh%state  = mesh%mstate+0.5d0*dt*mesh%hstate
      mesh%tstate = mesh%tstate+1.0d0/3.0d0*dt*mesh%hstate
    elseif(irk==3)then
      u  = mu+1.0d0*dt*hu
      tu = tu+1.0d0/3.0d0*dt*hu

      mesh%slip  = mesh%mslip+1.0d0*dt*mesh%sliprate
      mesh%tslip = mesh%tslip+1.0d0/3.0d0*dt*mesh%sliprate
      ! rate state
      mesh%state  = mesh%mstate+1.0d0*dt*mesh%hstate
      mesh%tstate = mesh%tstate+1.0d0/3.0d0*dt*mesh%hstate
    elseif(irk==4)then
      u  = tu+1.0d0/6.0d0*dt*hu

      mesh%slip  = mesh%tslip+1.0d0/6.0d0*dt*mesh%sliprate
      ! rate state
      mesh%state  = mesh%tstate+1.0d0/6.0d0*dt*mesh%hstate
    endif

  else if (timeIntegrationMethod == 1) then

    tu = rk4a(irk)*tu + dt*hu
    u = u + rk4b(irk)*tu

    mesh%tslip = rk4a(irk)*mesh%tslip + dt*mesh%sliprate
    mesh%slip = mesh%slip + rk4b(irk)*mesh%tslip
    ! rate state
    mesh%tstate = rk4a(irk)*mesh%tstate + dt*mesh%hstate
    mesh%state = mesh%state + rk4b(irk)*mesh%tstate
  else if (timeIntegrationMethod == 2) then

    tu = rk3a(irk)*tu + dt*hu
    u = u + rk3b(irk)*tu

    mesh%tslip = rk3a(irk)*mesh%tslip + dt*mesh%sliprate
    mesh%slip = mesh%slip + rk3b(irk)*mesh%tslip
    ! rate state
    mesh%tstate = rk3a(irk)*mesh%tstate + dt*mesh%hstate
    mesh%state = mesh%state + rk3b(irk)*mesh%tstate
  else if (timeIntegrationMethod == 3) then

    tu = rk2a(irk)*tu + dt*hu
    u = u + rk2b(irk)*tu

    mesh%tslip = rk2a(irk)*mesh%tslip + dt*mesh%sliprate
    mesh%slip = mesh%slip + rk2b(irk)*mesh%tslip
    ! rate state
    mesh%tstate = rk2a(irk)*mesh%tstate + dt*mesh%hstate
    mesh%state = mesh%state + rk2b(irk)*mesh%tstate

  end if

end subroutine

end module
