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

module mod_energy

  use mod_types, only : meshvar
  use mod_para,  only : Np, RKIND, CUSTOM_REAL, BC_FAULT, Nfaces, Nfp
  use mod_mpi,   only : sum_real
  use mod_eqns,  only : strain2stress

  implicit none

contains

subroutine cal_energy(mesh,u,Ev_sum,Es_sum)
  implicit none
  type(meshvar) :: mesh
  real(kind=RKIND),dimension(:,:) :: u
  integer :: ie,i,k1,k2
  real(kind=RKIND) :: Ev,Es
  real(kind=RKIND) :: Ev_sum,Es_sum
  real(kind=RKIND) :: Vx,Vy,Vz,strain(6),stress(6)
  real(kind=RKIND) :: rho,vp,vs,vol

  Ev = 0
  Es = 0
  do ie = 1,mesh%Nelem
    rho = mesh%rho(ie)
    vp = mesh%vp(ie)
    vs = mesh%vs(ie)
    vol = mesh%vol(ie)

    k1 = (ie-1)*Np+1
    k2 = ie*Np
    Vx = sum(u(k1:k2,1))/dble(Np)
    Vy = sum(u(k1:k2,2))/dble(Np)
    Vz = sum(u(k1:k2,3))/dble(Np)
    Ev = Ev + (Vx**2+Vy**2+Vz**2)/rho * vol

    do i = 1,6
      strain(i) = sum(u(k1:k2,3+i))/dble(Np)
    end do

    call strain2stress( &
        strain(1), &
        strain(2), &
        strain(3), &
        strain(4), &
        strain(5), &
        strain(6), &
        rho,vp,vs, &
        stress(1), &
        stress(2), &
        stress(3), &
        stress(4), &
        stress(5), &
        stress(6) )

    do i = 1,6
      Es = Es + strain(i) * stress(i) * vol
    end do

  end do ! elem

  call sum_real(Ev,Ev_sum,CUSTOM_REAL)
  call sum_real(Es,Es_sum,CUSTOM_REAL)

  !if (mesh%rank==0) &
  !print*,'rank=',mesh%rank,'Ev=',Ev_sum, 'Es=',Es_sum

end subroutine

subroutine cal_moment(mesh,Mmt,Mmtr)
  implicit none
  type(meshvar) :: mesh
  integer :: ie,ief,is
  real(kind=RKIND) :: Mmt,Mmtr
  real(kind=RKIND) :: Mmt_sum,Mmtr_sum
  real(kind=RKIND) :: rho,vs,area,mu,slip,rate

  Mmt = 0
  Mmtr = 0
  do ief = 1,mesh%nfault_elem
    ie = mesh%fault2wave(ief)
    rho = mesh%rho(ie)
    vs = mesh%vs(ie)
    mu = rho*vs*vs
    do is = 1,Nfaces
      if (mesh%bctype(is,ie) >= BC_FAULT ) then
        area = mesh%faultarea(is,ief)

        slip = sum(mesh%slip(:,is,ief))/dble(Nfp)
        rate = sum(mesh%sliprate(:,is,ief))/dble(Nfp)

        Mmt = Mmt + mu * area * slip
        Mmtr = Mmtr + mu * area * rate
      end if
    end do
  end do ! elem

  Mmt = Mmt * 0.5
  Mmtr = Mmtr * 0.5

  call sum_real(Mmt,Mmt_sum,CUSTOM_REAL)
  call sum_real(Mmtr,Mmtr_sum,CUSTOM_REAL)

end subroutine

end module
