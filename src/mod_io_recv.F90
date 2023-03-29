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

module mod_io_recv

  use mod_para,  only : RKIND,           &
                        Nfp, Nfaces,     &
                        NGLL,            &
                        BC_FREE,         &
                        BC_FAULT,        &
                        MAX_NUM_RECV,    &
                        NVAR_RECV,       &
                        data_dir
  use mod_types, only : meshvar
  use mod_interp, only: tri_interp
  use netcdf

  implicit none

  integer :: recv_ncid
  integer :: recv_varid(20)
  real(kind=rkind),dimension(:,:),allocatable :: val!(MAX_NUM_RECV,10)

  !real(kind=rkind),allocatable,dimension(:,:) :: PGVx,PGVy,PGVz,PGVh
  !real(kind=rkind),allocatable,dimension(:,:) :: Vx,Vy,Vz
  !real(kind=rkind),allocatable,dimension(:,:) :: Ux,Uy,Uz
  !real(kind=rkind),allocatable,dimension(:,:) :: exx,eyy,ezz,eyz,exz,exy

contains

subroutine check(status)
  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status)),' in file ',__FILE__,' line ',__LINE__
    stop "Stopped"
  end if
end subroutine

subroutine check2(status,msg)
  implicit none
  integer, intent ( in) :: status
  character(len=*) :: msg

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    print *, trim(msg)
    stop 110
  end if
end subroutine

subroutine recv_io_init(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: myrank
  !integer :: i, j, ie, ief, is, ierr
  integer :: ierr
  character(len=128) :: filename
  !real(kind=rkind),allocatable,dimension(:,:) :: x,y,z
  !real(kind=rkind),allocatable,dimension(:,:,:) :: vars
  integer :: dimid(20),varid(20)

  if (mesh%nrecv==0) return

  myrank = mesh%rank

  write(filename,'(a,a,i6.6,a)') trim(data_dir),'/recv_mpi',myrank,'.nc'

  ! start to create
  ierr = nf90_create(filename, NF90_CLOBBER, recv_ncid)
  call check2(ierr,'nf90_create recv')

  ! define dimensions
  ierr = nf90_def_dim(recv_ncid, "Nrecv", mesh%nrecv, dimid(1))
  call check2(ierr,'def_dim nrecv')
  ierr = nf90_def_dim(recv_ncid, "Nt", nf90_unlimited, dimid(2))
  call check2(ierr,'def_dim Nt')
  ierr = nf90_def_dim(recv_ncid, "three", 3, dimid(3))
  call check2(ierr,'def_dim 3')
  ierr = nf90_def_dim(recv_ncid, "Nvar", NVAR_RECV, dimid(4))
  call check2(ierr,'def_dim Nvar')

  ! define variables
  ierr = nf90_def_var(recv_ncid, "id", NF90_INT, dimid(1), varid(1))
  call check2(ierr,'def_var recv_id')
  ierr = nf90_def_var(recv_ncid, "bctype", NF90_INT, dimid(1), varid(2))
  call check2(ierr,'def_var bctype')
  ierr = nf90_def_var(recv_ncid, "coord", NF90_DOUBLE, (/dimid(3),dimid(1)/), varid(3))
  call check2(ierr,'def_var coord')
  ierr = nf90_def_var(recv_ncid, "normal", NF90_DOUBLE, (/dimid(3),dimid(1)/), varid(4))
  call check2(ierr,'def_var normal')

  ierr = nf90_def_var(recv_ncid, "time", NF90_DOUBLE, dimid(2), recv_varid(1))
  call check2(ierr,'def_var time')
  ierr = nf90_def_var(recv_ncid, "var", NF90_FLOAT, (/dimid(1),dimid(4),dimid(2)/), recv_varid(2))
  call check2(ierr,'def_var vars')

  ! end of define
  ierr = nf90_enddef(recv_ncid)
  call check2(ierr,'enddef')

  ! put variables
  ierr = nf90_put_var(recv_ncid, varid(1), mesh%recv_id)
  call check2(ierr,'put_var recv_id')
  ierr = nf90_put_var(recv_ncid, varid(2), mesh%recv_bctype)
  call check2(ierr,'put_var recv_bctype')
  ierr = nf90_put_var(recv_ncid, varid(3), mesh%recv_coord)
  call check2(ierr,'put_var recv_coord')
  ierr = nf90_put_var(recv_ncid, varid(4), mesh%recv_normal)
  call check2(ierr,'put_var recv_normal')

  allocate(val(mesh%nrecv,NVAR_RECV))
end subroutine

subroutine recv_io_save(mesh,u,it)
  implicit none
  type(meshvar) :: mesh
  real(kind=rkind) :: u(:,:)
  integer :: i, j, ie, ief, is, it, ierr, n
  !integer,dimension(3) :: start,cnt,stride
  !real(kind=rkind) :: val(MAX_NUM_RECV,10)
  real(kind=rkind),dimension(3) :: v1,v2,v3,p
  real(kind=rkind) :: c1,c2,c3
  if (mesh%nrecv==0) return

  ierr = nf90_put_var(recv_ncid,recv_varid(1),(/mesh%current_time/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var time')

  do ie = 1,mesh%nelem
    do is = 1,Nfaces
      do n = 1,mesh%nrecv
        if ( mesh%recv_bctype(n) == BC_FREE .and. &
            ie == mesh%recv_elem(n) .and. is == mesh%recv_face(n)) then
          mesh%recv_buffer(:,n,1) = u(mesh%vmapM(:,is,ie),1)/mesh%rho(ie)
          mesh%recv_buffer(:,n,2) = u(mesh%vmapM(:,is,ie),2)/mesh%rho(ie)
          mesh%recv_buffer(:,n,3) = u(mesh%vmapM(:,is,ie),3)/mesh%rho(ie)
        end if
        if ( mesh%recv_bctype(n) >= BC_FAULT .and. &
            ie == mesh%recv_elem(n) .and. is == mesh%recv_face(n)) then
          ief = mesh%wave2fault(ie)
          mesh%recv_buffer(:,n,1)  = mesh%sliprate1(:,is,ief)
          mesh%recv_buffer(:,n,2)  = mesh%sliprate2(:,is,ief)
          mesh%recv_buffer(:,n,3)  = mesh%stress1  (:,is,ief)
          mesh%recv_buffer(:,n,4)  = mesh%stress2  (:,is,ief)
          mesh%recv_buffer(:,n,5)  = mesh%sigma    (:,is,ief)
          mesh%recv_buffer(:,n,6)  = mesh%slip1    (:,is,ief)
          mesh%recv_buffer(:,n,7)  = mesh%slip2    (:,is,ief)
          mesh%recv_buffer(:,n,8)  = mesh%state    (:,is,ief)
          mesh%recv_buffer(:,n,9)  = mesh%TP_T     (:,is,ief)
          mesh%recv_buffer(:,n,10) = mesh%TP_P     (:,is,ief)
          mesh%recv_buffer(:,n,11) = u(mesh%vmapM(:,is,ie),1)/mesh%rho(ie)
          mesh%recv_buffer(:,n,12) = u(mesh%vmapM(:,is,ie),2)/mesh%rho(ie)
          mesh%recv_buffer(:,n,13) = u(mesh%vmapM(:,is,ie),3)/mesh%rho(ie)
        end if
      end do
    end do
  enddo

  do j = 1,mesh%nrecv
    do i = 1,NVAR_RECV
      v1=(/ &
          mesh%vx(mesh%vmapM(1,mesh%recv_face(j),mesh%recv_elem(j))), &
          mesh%vy(mesh%vmapM(1,mesh%recv_face(j),mesh%recv_elem(j))), &
          mesh%vz(mesh%vmapM(1,mesh%recv_face(j),mesh%recv_elem(j))) /)
      v2=(/ &
          mesh%vx(mesh%vmapM(NGLL,mesh%recv_face(j),mesh%recv_elem(j))), &
          mesh%vy(mesh%vmapM(NGLL,mesh%recv_face(j),mesh%recv_elem(j))), &
          mesh%vz(mesh%vmapM(NGLL,mesh%recv_face(j),mesh%recv_elem(j))) /)
      v3=(/ &
          mesh%vx(mesh%vmapM(Nfp,mesh%recv_face(j),mesh%recv_elem(j))), &
          mesh%vy(mesh%vmapM(Nfp,mesh%recv_face(j),mesh%recv_elem(j))), &
          mesh%vz(mesh%vmapM(Nfp,mesh%recv_face(j),mesh%recv_elem(j))) /)

      c1 = mesh%recv_buffer(1   ,j,i)
      c2 = mesh%recv_buffer(NGLL,j,i)
      c3 = mesh%recv_buffer(Nfp ,j,i)

      p = (/ &
          mesh%recv_coord(1,j), &
          mesh%recv_coord(2,j), &
          mesh%recv_coord(3,j) /)

      val(j,i) = tri_interp(v1,v2,v3,c1,c2,c3,p)
      !val(j,i) = interp_dist( &
      !mesh%vx(mesh%vmapM(:,mesh%recv_face(j),mesh%recv_elem(j))), &
      !mesh%vy(mesh%vmapM(:,mesh%recv_face(j),mesh%recv_elem(j))), &
      !mesh%vz(mesh%vmapM(:,mesh%recv_face(j),mesh%recv_elem(j))), &
      !mesh%recv_buffer(:,j,i), &
      !Nfp, &
      !mesh%recv_coord(1,j), &
      !mesh%recv_coord(2,j), &
      !mesh%recv_coord(3,j) )
    end do
  end do

  !val(1) = mesh%recv_buffer(1,1,1)

  ierr = nf90_put_var(recv_ncid,recv_varid(2), &
      (/sngl(val(1:mesh%nrecv,1:NVAR_RECV))/),(/1,1,it/),(/mesh%nrecv,NVAR_RECV,1/),(/1,1/))
  call check2(ierr,'put_var vars')

  ierr = nf90_sync(recv_ncid)

end subroutine

subroutine recv_io_end(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: ierr
  if (mesh%nrecv==0) return

  ierr = nf90_close(recv_ncid)
  call check2(ierr,'nf90_close recv')

  if(allocated(val)) deallocate(val)

end subroutine

end module
