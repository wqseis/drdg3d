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

module mod_recv

  use mod_para, only : RKIND,         &
                       MAX_NUM_RECV,  &
                       BC_FREE,       &
                       BC_FAULT,      &
                       Nfp,           &
                       Nfaces
  use mod_math, only : cross
  use mod_mesh, only : meshvar

  implicit none

contains

subroutine locate_recvs(mesh)
  implicit none

  type(meshvar) :: mesh
  integer :: ir,ie,is,i,n
  !character(len=80) :: filename,namestr,tmpstr
  character(len=80) :: tmpstr
  integer :: ierr!,ncid,varid
  integer :: nr,myrank
  integer :: fid

  real(kind=rkind),dimension(3) :: v1,v2,v3,v4,p
  real(kind=rkind),dimension(3,MAX_NUM_RECV) :: recv_coord
  real(kind=rkind),dimension(3,MAX_NUM_RECV) :: recv_coord_loc

  integer,dimension(MAX_NUM_RECV) :: recv_bctype, recv_elem, recv_face

  myrank = mesh%rank

  !recv_coord(:,1) = (/1e-16,0.0,-7.5/)
  !recv_coord(:,2) = (/1e-16,3.0,-7.5/)
  !recv_coord(:,3) = (/1e-16,4.5,-7.5/)
  !recv_coord(:,4) = (/1e-16,7.5,-7.5/)
  !recv_coord(:,5) = (/1e-16,9.0,-7.5/)
  !recv_coord(:,6) = (/1e-16,12.,-7.5/)
  !recv_coord(:,7) = (/1e-16,0.0,-3.0/)
  !recv_coord(:,8) = (/1e-16,0.0,-1.5/)
  !recv_coord(:,9) = (/1e-16,0.0,-12.0/)
  fid = 100
  n = -100
  open(fid,file='stations.txt',iostat=ierr,status='old')
  if(ierr .ne. 0) then
    if(myrank==0) then
      print*,'Warning! I cannot open stations.txt. Set num of recv = 0'
      n = 0
    end if
  else
    read(fid,'(a)') tmpstr
    read(fid,*) n
    !print*,tmpstr
    !print*,n
    do i = 1,n
      read(fid,*) recv_coord(1,i),recv_coord(2,i),recv_coord(3,i),recv_bctype(i)
      !print*,recv_coord(1,i),recv_coord(2,i),recv_coord(3,i),recv_bctype(i)
    end do
  end if
  close(fid)

  nr = 0

  do ir = 1,n!mesh%nrecv
    do ie = 1,mesh%nelem
      v1=mesh%coord(:,mesh%elem(1,ie))
      v2=mesh%coord(:,mesh%elem(2,ie))
      v3=mesh%coord(:,mesh%elem(3,ie))
      v4=mesh%coord(:,mesh%elem(4,ie))
      p = recv_coord(:,ir)
      if(PointInTet(v1,v2,v3,v4,p)) then
        nr=nr+1
        do is = 1,Nfaces
          if(mesh%bctype(is,ie) >= BC_FAULT) then
            !print*,'rank=',myrank,'recv elem=',ie,'face=',is, 'coord=', &
            !sngl(mesh%vx(mesh%vmapM(:,is,ie))),&
            !sngl(mesh%vy(mesh%vmapM(:,is,ie))),&
            !sngl(mesh%vz(mesh%vmapM(:,is,ie)))
            recv_elem(nr) = ie
            recv_face(nr) = is
            recv_coord_loc(:,nr) = p
            recv_bctype(nr) = BC_FAULT
          end if
          if(mesh%bctype(is,ie) == BC_FREE) then
            !print*,'rank=',myrank,'recv elem=',ie,'face=',is, 'coord=', &
            !sngl(mesh%vx(mesh%vmapM(:,is,ie))),&
            !sngl(mesh%vy(mesh%vmapM(:,is,ie))),&
            !sngl(mesh%vz(mesh%vmapM(:,is,ie)))
            recv_elem(nr) = ie
            recv_face(nr) = is
            recv_coord_loc(:,nr) = p
            recv_bctype(nr) = BC_FREE
          end if
        end do
      end if
    end do
  end do

  mesh%nrecv = nr
  if (mesh%nrecv>0) then
    allocate(mesh%recv_bctype(1:mesh%nrecv))
    allocate(mesh%recv_elem(1:mesh%nrecv))
    allocate(mesh%recv_face(1:mesh%nrecv))
    allocate(mesh%recv_coord(1:3,1:mesh%nrecv))
    allocate(mesh%recv_buffer(Nfp,mesh%nrecv,10))
    mesh%recv_bctype = recv_bctype(1:mesh%nrecv)
    mesh%recv_elem = recv_elem(1:mesh%nrecv)
    mesh%recv_face = recv_face(1:mesh%nrecv)
    mesh%recv_coord = recv_coord_loc(:,1:mesh%nrecv)
    mesh%recv_buffer = 0.
  end if



  print*,'rank = ',myrank,'nrecv = ',mesh%nrecv


end subroutine

! flag = PointInTet(v1,v2,v3,v4,p)
! https://people.sc.fsu.edu/~jburkardt/presentations/cg_lab_barycentric_tetrahedrons.pdf
function PointInTet(v1,v2,v3,v4,p)

  implicit none

  logical PointInTet

  real(kind=RKIND),dimension(3) :: v1,v2,v3,v4,p
  real(kind=RKIND) :: c1,c2,c3,c4

  PointInTet = .False.

  call barycentricPointInTet(v1,v2,v3,v4,p,c1,c2,c3,c4)

  if(c1>=0.0 .and. c2>=0.0 .and. c3>=0.0 .and. c4>=0.0) then
    PointInTet = .True.
  end if

end function

subroutine barycentricPointInTet(v1,v2,v3,v4,p,c1,c2,c3,c4)
  implicit none

  real(kind=RKIND),dimension(3) :: v1,v2,v3,v4,p
  real(kind=RKIND) :: c1,c2,c3,c4

  c1 = SignedDistance(v2,v3,v4,p) / &
       SignedDistance(v2,v3,v4,v1)
  c2 = SignedDistance(v1,v3,v4,p) / &
       SignedDistance(v1,v3,v4,v2)
  c3 = SignedDistance(v1,v2,v4,p) / &
       SignedDistance(v1,v2,v4,v3)
  c4 = SignedDistance(v1,v2,v3,p) / &
       SignedDistance(v1,v2,v3,v4)

end

function SignedDistance(v1,v2,v3,p)

  implicit none

  real(kind=rkind) :: SignedDistance
  real(kind=RKIND),dimension(3) :: v1,v2,v3,p,n

  n = cross(v2-v1,v3-v1)
  SignedDistance = dot_product(n, p-v1)

end function

end module
