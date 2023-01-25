!--------------------------------------------------------------------------
!   Copyright 2012-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2015-2017 Andre Lamert (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of NEXD 3D.
!
!   NEXD 3D is free software: you can redistribute it and/or modify it
!   under the terms of the GNU General Public License as published by the
!   Free Software Foundation, either version 3 of the License, or (at your
!   option) any later version.
!
!   NEXD 3D is distributed in the hope that it will be useful, but WITHOUT
!   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!   FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License v3.0
!   along with NEXD 3D. If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------
module mod_vtk
! module to create files in VTK format which can be read by paraview
  use mod_para, only : RKIND, CUSTOM_REAL
  implicit none

contains

subroutine writeVtkTetraMeshIntdata(filename, elem, coord, data)
  implicit none
  character(len=*) :: filename
  real(kind=CUSTOM_REAL), dimension(:,:) :: coord
  integer(kind=CUSTOM_REAL), dimension(:) :: data
  integer, dimension(:,:) :: elem
  integer :: nelem, ncoord
  integer :: i

  nelem=size(elem(1,:))
  ncoord=size(coord(1,:))

  open(19,file=trim(filename),status='unknown')
  write(19,'(a)') '# vtk DataFile Version 3.1'
  write(19,'(a)') 'Mesh with Data'
  write(19,'(a)') 'ASCII'
  write(19,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(19, '(a,i12,a)') 'POINTS ', ncoord, ' float'
  do i=1,ncoord
    write(19,'(3e18.6)') coord(1,i),coord(2,i),coord(3,i)
  enddo
  write(19,*) ""

  write(19,'(a,i12,i12)') "CELLS ",nelem,nelem*5 
  do i=1,nelem
    write(19,'(5i12)') 4 ,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1,elem(4,i)-1
  end do
  write(19,*) ""

  write(19,'(a,i12)') "CELL_TYPES ",nelem
  write(19,*) (10,i=1,nelem)
  write(19,*) ""

  write(19,'(a,i12)') "CELL_DATA ",nelem
  write(19,'(a)') "SCALARS data integer"
  write(19,'(a)') "LOOKUP_TABLE default"
  do i = 1,nelem
    write(19,*) data(i)
  enddo
  write(19,*) ""
  close(19)

end subroutine writeVtkTetraMeshIntdata

subroutine writeVtkTetraMeshRealdata(filename, elem, coord, data)
  implicit none
  character(len=*) :: filename
  real(kind=CUSTOM_REAL), dimension(:,:) :: coord
  real(kind=CUSTOM_REAL), dimension(:) :: data
  integer, dimension(:,:) :: elem
  integer :: nelem, ncoord
  integer :: i

  nelem=size(elem(1,:))
  ncoord=size(coord(1,:))

  open(19,file=trim(filename),status='unknown')
  write(19,'(a)') '# vtk DataFile Version 3.1'
  write(19,'(a)') 'Mesh with Data'
  write(19,'(a)') 'ASCII'
  write(19,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(19, '(a,i12,a)') 'POINTS ', ncoord, ' float'
  do i=1,ncoord
    write(19,'(3e18.6)') coord(1,i),coord(2,i),coord(3,i)
  enddo
  write(19,*) ""

  write(19,'(a,i12,i12)') "CELLS ",nelem,nelem*5 
  do i=1,nelem
    write(19,'(5i12)') 4 ,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1,elem(4,i)-1
  end do
  write(19,*) ""

  write(19,'(a,i12)') "CELL_TYPES ",nelem
  write(19,*) (10,i=1,nelem)
  write(19,*) ""

  write(19,'(a,i12)') "CELL_DATA ",nelem
  write(19,'(a)') "SCALARS data float"
  write(19,'(a)') "LOOKUP_TABLE default"
  do i = 1,nelem
    write(19,*) data(i)
  enddo
  write(19,*) ""
  close(19)

end subroutine writeVtkTetraMeshRealdata

subroutine writeVtkTetraNodeRealdata(filename, elem, coord, data)
  implicit none
  character(len=*) :: filename
  real(kind=CUSTOM_REAL), dimension(:,:) :: coord
  real(kind=CUSTOM_REAL), dimension(:) :: data
  integer, dimension(:,:) :: elem
  integer :: nelem, ncoord
  integer :: i

  nelem=size(elem(1,:))
  ncoord=size(coord(1,:))

  open(19,file=trim(filename),status='unknown')
  write(19,'(a)') '# vtk DataFile Version 3.1'
  write(19,'(a)') 'Mesh with Data'
  write(19,'(a)') 'ASCII'
  write(19,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(19, '(a,i12,a)') 'POINTS ', ncoord, ' float'
  do i=1,ncoord
    write(19,'(3e18.6)') coord(1,i),coord(2,i),coord(3,i)
  enddo
  write(19,*) ""

  write(19,'(a,i12,i12)') "CELLS ",nelem,nelem*5 
  do i=1,nelem
    write(19,'(5i12)') 4 ,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1,elem(4,i)-1
  end do
  write(19,*) ""

  write(19,'(a,i12)') "CELL_TYPES ",nelem
  write(19,*) (10,i=1,nelem)
  write(19,*) ""

  write(19,'(a,i12)') "CELL_DATA ",nelem
  write(19,'(a)') "SCALARS data float"
  write(19,'(a)') "LOOKUP_TABLE default"
  do i = 1,nelem
    write(19,*) data(i)
  enddo
  write(19,*) ""
  close(19)

end subroutine writeVtkTetraNodeRealdata

subroutine writeVtkNodes(filename, x,y,z)
  implicit none
  character(len=*) :: filename
  real(kind=CUSTOM_REAL), dimension(:) ::x,y,z
  integer :: n
  integer :: i

  n=size(x)

  open(19,file=trim(filename),status='unknown')
  write(19,'(a)') '# vtk DataFile Version 3.1'
  write(19,'(a)') 'Mesh with Data'
  write(19,'(a)') 'ASCII'
  write(19,'(a)') 'DATASET POLYDATA'
  write(19, '(a,i12,a)') 'POINTS ', n, ' float'
  do i=1,n
    write(19,'(3e18.6)') x(i),y(i),z(i)
  enddo
  write(19,*) ""
  write(19,'(a,i12,i12)') "VERTICES ",n,n*4 
  do i=1,n
    write(19,'(5i12)') 3 ,i-1,i-1,i-1
  end do
  write(19,*) ""

  close(19)

end subroutine writeVtkNodes

subroutine writeVtkNodesRealData(filename, x,y,z,data)
  implicit none
  character(len=*) :: filename
  real(kind=CUSTOM_REAL), dimension(:) ::x,y,z,data
  integer :: n
  integer :: i

  n=size(x)

  open(19,file=trim(filename),status='unknown')
  write(19,'(a)') '# vtk DataFile Version 3.1'
  write(19,'(a)') 'Mesh with Data'
  write(19,'(a)') 'ASCII'
  write(19,'(a)') 'DATASET POLYDATA'
  write(19, '(a,i12,a)') 'POINTS ', n, ' float'
  do i=1,n
    write(19,'(3e18.6)') x(i),y(i),z(i)
  enddo
  write(19,*) ""
  write(19,'(a,i12,i12)') "VERTICES ",n,n*4 
  do i=1,n
    write(19,'(5i12)') 3 ,i-1,i-1,i-1
  end do
  write(19,'(a,i12)') "POINT_DATA ",n
  write(19,'(a)') "SCALARS data float"
  write(19,'(a)') "LOOKUP_TABLE data"
  do i = 1,n
    write(19,*) data(i)
  enddo
  write(19,*) ""
  close(19)

end subroutine writeVtkNodesRealData

subroutine writeVtkNodesRealDataBin(filename, x,y,z,data)
  implicit none
  character(len=*) :: filename
  real(kind=CUSTOM_REAL), dimension(:) ::x,y,z,data
  integer :: n
  integer :: i

  n=size(x)

  open(19,file=trim(filename),status='unknown',form='UNFORMATTED')
  write(19) '# vtk DataFile Version 3.1'
  write(19) 'Mesh with Data'
  write(19) 'ASCII'
  write(19) 'DATASET POLYDATA'
  write(19) 'POINTS ', n, ' float'
  do i=1,n
    write(19) x(i),y(i),z(i)
  enddo
  write(19) ""
  write(19) "VERTICES ",n,n*4 
  do i=1,n
    write(19) 3 ,i-1,i-1,i-1
  end do
  write(19) "POINT_DATA ",n
  write(19) "SCALARS data float"
  write(19) "LOOKUP_TABLE data"
  do i = 1,n
    write(19) data(i)
  enddo
  write(19) ""
  close(19)

end subroutine writeVtkNodesRealDataBin

end module
