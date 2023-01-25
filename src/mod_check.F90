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

module mod_check

  use mod_types,  only : meshvar

  implicit none

contains

subroutine check_geometry(mesh)
  implicit none
  type(meshvar) :: mesh

  print*,'Fmask=',mesh%Fmask
  print*,'LIFT=',sngl(transpose(mesh%LIFT))
  print*,'rx = ',sngl(minval(mesh%rx)),'~',sngl(maxval(mesh%rx))
  print*,'ry = ',sngl(minval(mesh%ry)),'~',sngl(maxval(mesh%ry))
  print*,'rz = ',sngl(minval(mesh%rz)),'~',sngl(maxval(mesh%rz))
  print*,'sx = ',sngl(minval(mesh%sx)),'~',sngl(maxval(mesh%sx))
  print*,'sy = ',sngl(minval(mesh%sy)),'~',sngl(maxval(mesh%sy))
  print*,'sz = ',sngl(minval(mesh%sz)),'~',sngl(maxval(mesh%sz))
  print*,'tx = ',sngl(minval(mesh%tx)),'~',sngl(maxval(mesh%tx))
  print*,'ty = ',sngl(minval(mesh%ty)),'~',sngl(maxval(mesh%ty))
  print*,'tz = ',sngl(minval(mesh%tz)),'~',sngl(maxval(mesh%tz))
  print*,'J = ',sngl(minval(mesh%jac)),'~',sngl(maxval(mesh%jac))
  print*,'sJ = ',sngl(minval(mesh%sj)),'~',sngl(maxval(mesh%sj))
  print*,'Fscale = ',sngl(minval(mesh%Fscale)),'~',sngl(maxval(mesh%Fscale))
  print*,'nx = ',sngl(minval(mesh%nx)),'~',sngl(maxval(mesh%nx))
  print*,'ny = ',sngl(minval(mesh%ny)),'~',sngl(maxval(mesh%ny))
  print*,'nz = ',sngl(minval(mesh%nz)),'~',sngl(maxval(mesh%nz))
  print*,'mx = ',sngl(minval(mesh%mx)),'~',sngl(maxval(mesh%mx))
  print*,'my = ',sngl(minval(mesh%my)),'~',sngl(maxval(mesh%my))
  print*,'mz = ',sngl(minval(mesh%mz)),'~',sngl(maxval(mesh%mz))
  print*,'lx = ',sngl(minval(mesh%lx)),'~',sngl(maxval(mesh%lx))
  print*,'ly = ',sngl(minval(mesh%ly)),'~',sngl(maxval(mesh%ly))
  print*,'lz = ',sngl(minval(mesh%lz)),'~',sngl(maxval(mesh%lz))
end subroutine

end module
