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

module mod_damp
  use mod_para, only : RKIND, CUSTOM_REAL, Np
  use mod_mesh, only : meshvar
  use mod_mpi,  only : minval_real_all, &
                       maxval_real_all
  use mpi
  implicit none

contains

subroutine cal_coord_range(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: rank
  real(kind=RKIND) :: xmin,ymin,zmin
  real(kind=RKIND) :: xmax,ymax,zmax

  rank = mesh%rank

  xmin = minval(mesh%vx)
  ymin = minval(mesh%vy)
  zmin = minval(mesh%vz)
  xmax = maxval(mesh%vx)
  ymax = maxval(mesh%vy)
  zmax = maxval(mesh%vz)

  call minval_real_all(xmin,mesh%xmin,CUSTOM_REAL)
  call minval_real_all(ymin,mesh%ymin,CUSTOM_REAL)
  call minval_real_all(zmin,mesh%zmin,CUSTOM_REAL)
  call maxval_real_all(xmax,mesh%xmax,CUSTOM_REAL)
  call maxval_real_all(ymax,mesh%ymax,CUSTOM_REAL)
  call maxval_real_all(zmax,mesh%zmax,CUSTOM_REAL)

  if(rank==0) print*,'x = ',sngl(mesh%xmin),'~',sngl(mesh%xmax)
  if(rank==0) print*,'y = ',sngl(mesh%ymin),'~',sngl(mesh%ymax)
  if(rank==0) print*,'z = ',sngl(mesh%zmin),'~',sngl(mesh%zmax)

end subroutine

subroutine init_damp(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: i,j,ie
  integer :: nelem,rank
  real(kind=rkind) :: x1,x2,y1,y2,z1,z2,dampLen,dpx,dpy,dpz
  !real(kind=rkind) :: r,r1,x,y,z,alpha,coef
  real(kind=rkind) :: x,y,z,alpha,coef
  !real(kind=rkind) :: axisdampLenimits(6)
  !character(len=80) :: filename

  rank = mesh%rank


  coef = 0.92d0
  dampLen = 20d0 ! km
  alpha = sqrt(-log(coef))
  if(rank==0) print*,'damp coeff = ', coef
  if(rank==0) print*,'damp Len   = ', dampLen
  if(rank==0) print*,'damp alpha = ', alpha

  x1 = mesh%xmin
  x2 = mesh%xmax
  y1 = mesh%ymin
  y2 = mesh%ymax
  z1 = mesh%zmin
  z2 = mesh%zmax

  !x1 = -50
  !x2 = 50
  !y1 = -50
  !y2 = 50
  !z1 = -50
  !z2 = 0

  dpx = 1d0
  dpy = 1d0
  dpz = 1d0

  nelem = mesh%nelem
  
  allocate(mesh%damp(Np*Nelem))

  do ie = 1,nelem
    do i = 1,Np
      j = i+(ie-1)*Np
      x = mesh%vx(j)
      y = mesh%vy(j)
      z = mesh%vz(j)

      if (x < x1+dampLen) then
          !dpx = dexp(-0.1*(x-(x1+dampLen))**2/dampLen**2)
          dpx = dexp(-(alpha*(x-(x1+dampLen))/dampLen)**2)
          !print*,'damp = ', x-(x1+dampLen),dpx
      else if (x >= x1+dampLen .and. x <= x2-dampLen) then
          dpx = 1d0
      else
          !dpx = dexp(-0.1*(x-(x2-dampLen))**2/dampLen**2)
          dpx = dexp(-(alpha*(x-(x2-dampLen))/dampLen)**2)
      end if
  
      if (y < y1+dampLen) then
          !dpy = dexp(-0.1*(y-(y1+dampLen))**2/dampLen**2)
          dpy = dexp(-(alpha*(y-(y1+dampLen))/dampLen)**2)
      else if (y >= y1+dampLen .and. y <= y2-dampLen) then
          dpy = 1d0
      else
          !dpy = dexp(-0.1*(y-(y2-dampLen))**2/dampLen**2)
          dpy = dexp(-(alpha*(y-(y2-dampLen))/dampLen)**2)
      end if

      if (z < z1+dampLen) then
          !dpz = dexp(-0.1*(z-(z1+dampLen))**2/dampLen**2)
          dpz = dexp(-(alpha*(z-(z1+dampLen))/dampLen)**2)
      else! if (z >= z1+dampLen .and. z <= z2-dampLen) then
          dpz = 1d0
      !else
      !    dpz = dexp(-0.1*(z-(z2-dampLen))**2/dampLen**2)
      end if

      mesh%damp(j) = min(min(dpx,dpy),dpz)

    end do
  end do

end subroutine

end module
