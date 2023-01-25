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
                        BC_FREE,         &
                        BC_FAULT,        &
                        MAX_NUM_RECV,    &
                        data_dir
  use mod_types, only : meshvar
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
  !print*,filename
  !allocate(x(Nfp,mesh%nfree_face))
  !allocate(y(Nfp,mesh%nfree_face))
  !allocate(z(Nfp,mesh%nfree_face))

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
  ierr = nf90_def_dim(recv_ncid, "ten", 10, dimid(4))
  call check2(ierr,'def_dim 10')

  ! define variables
  ierr = nf90_def_var(recv_ncid, "coord" , NF90_DOUBLE, (/dimid(3),dimid(1)/), varid(1))
  call check2(ierr,'def_var coord')
  ierr = nf90_def_var(recv_ncid, "bctype" , NF90_INT, dimid(1), varid(2))
  call check2(ierr,'def_var bctype')

  ierr = nf90_def_var(recv_ncid, "time", NF90_DOUBLE, dimid(2), recv_varid(1))
  call check2(ierr,'def_var time')
  ierr = nf90_def_var(recv_ncid, "var", NF90_FLOAT, (/dimid(1),dimid(4),dimid(2)/), recv_varid(2))
  call check2(ierr,'def_var vars')
  !ierr = nf90_def_var(recv_ncid, "sliprate", NF90_FLOAT, dimid(1:2), recv_varid(2))
  !call check2(ierr,'def_var sliprate')
  !ierr = nf90_def_var(recv_ncid, "stress", NF90_FLOAT, dimid(1:2), recv_varid(3))
  !call check2(ierr,'def_var stress')
  !ierr = nf90_def_var(recv_ncid, "sigma", NF90_FLOAT, dimid(1:2), recv_varid(4))
  !call check2(ierr,'def_var sigma')
  !ierr = nf90_def_var(recv_ncid, "slip", NF90_FLOAT, dimid(1:2), recv_varid(5))
  !call check2(ierr,'def_var slip')

  ! end of define
  ierr = nf90_enddef(recv_ncid)
  call check2(ierr,'enddef')

  ! put variables
  ierr = nf90_put_var(recv_ncid, varid(1), mesh%recv_coord)
  call check2(ierr,'put_var recv_coord')
  ierr = nf90_put_var(recv_ncid, varid(2), mesh%recv_bctype)
  call check2(ierr,'put_var recv_bctype')
  !ierr = nf90_put_var(recv_ncid, varid(2), y)
  !call check2(ierr,'put_var y')
  !ierr = nf90_put_var(recv_ncid, varid(3), z)
  !call check2(ierr,'put_var z')

  !!!ierr = nf90_close(recv_ncid)
  !deallocate(x,y,z)
  !deallocate(vars)
  allocate(val(mesh%nrecv,10))
end subroutine

subroutine recv_io_save(mesh,u,it)
  implicit none
  type(meshvar) :: mesh
  real(kind=rkind) :: u(:,:)
  integer :: i, j, ie, ief, is, it, ierr, n
  !integer,dimension(3) :: start,cnt,stride
  !real(kind=rkind) :: val(MAX_NUM_RECV,10)
  if (mesh%nrecv==0) return

  !i = 0
  !do ie = 1,mesh%nelem
  !  do is = 1,Nfaces
  !    if ( mesh%bctype(is,ie) == BC_FREE) then
  !      i = i + 1
  !      Vx(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),1)/mesh%rho(ie) )
  !      Vy(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),2)/mesh%rho(ie) )
  !      Vz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),3)/mesh%rho(ie) )
  !      exx(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),4))
  !      eyy(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),5))
  !      ezz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),6))
  !      eyz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),4))
  !      exz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),5))
  !      exy(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),6))
  !      Ux(1:Nfp,i) = sngl(displ(mesh%vmapM(1:Nfp,is,ie),1)/mesh%rho(ie) )
  !      Uy(1:Nfp,i) = sngl(displ(mesh%vmapM(1:Nfp,is,ie),2)/mesh%rho(ie) )
  !      Uz(1:Nfp,i) = sngl(displ(mesh%vmapM(1:Nfp,is,ie),3)/mesh%rho(ie) )
  !    endif
  !  enddo
  !enddo
  ! PGV
  !do i = 1,Nfp
  !  do j = 1,mesh%nfree_face
  !    PGVx(i,j) = max(PGVx(i,j),abs(Vx(i,j)))
  !    PGVy(i,j) = max(PGVy(i,j),abs(Vy(i,j)))
  !    PGVz(i,j) = max(PGVz(i,j),abs(Vz(i,j)))
  !  enddo
  !end do

  !start=(/1,1,it+0/); cnt=(/Nfp,mesh%nfree_face,1/); stride=(/1,1,1/)
  !ierr = nf90_put_var(recv_ncid,recv_varid(1),Vx,start,cnt,stride)
  !call check2(ierr,'put_var Vx')
  !ierr = nf90_put_var(recv_ncid,recv_varid(2),Vy,start,cnt,stride)
  !call check2(ierr,'put_var Vy')
  !ierr = nf90_put_var(recv_ncid,recv_varid(3),Vz,start,cnt,stride)
  !call check2(ierr,'put_var Vz')


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
        if ( mesh%recv_bctype(n) == BC_FAULT .and. &
            ie == mesh%recv_elem(n) .and. is == mesh%recv_face(n)) then
          ief = mesh%wave2fault(ie)
          mesh%recv_buffer(:,n,1) = mesh%sliprate(:,is,ief)
          mesh%recv_buffer(:,n,2) = mesh%stress  (:,is,ief)
          mesh%recv_buffer(:,n,3) = mesh%sigma   (:,is,ief)
          mesh%recv_buffer(:,n,4) = mesh%slip    (:,is,ief)
        end if
      end do
    end do
  enddo

  do j = 1,mesh%nrecv
    do i = 1,4
      val(j,i) = interp_dist( &
      mesh%vx(mesh%vmapM(:,mesh%recv_face(j),mesh%recv_elem(j))), &
      mesh%vy(mesh%vmapM(:,mesh%recv_face(j),mesh%recv_elem(j))), &
      mesh%vz(mesh%vmapM(:,mesh%recv_face(j),mesh%recv_elem(j))), &
      mesh%recv_buffer(:,j,i), &
      Nfp, &
      mesh%recv_coord(1,j), &
      mesh%recv_coord(2,j), &
      mesh%recv_coord(3,j) )
    end do
  end do

  !val(1) = mesh%recv_buffer(1,1,1)

  ierr = nf90_put_var(recv_ncid,recv_varid(2),(/sngl(val(1:mesh%nrecv,1:10))/),(/1,1,it/),(/mesh%nrecv,10,1/),(/1,1/))
  call check2(ierr,'put_var vars')
  !ierr = nf90_put_var(recv_ncid,recv_varid(2),(/sngl(val(1:mesh%nrecv,1))/),(/1,it/),(/mesh%nrecv,1/),(/1,1/))
  !call check2(ierr,'put_var sliprate')
  !ierr = nf90_put_var(recv_ncid,recv_varid(3),(/sngl(val(1:mesh%nrecv,2))/),(/1,it/),(/mesh%nrecv,1/),(/1,1/))
  !call check2(ierr,'put_var stress')
  !ierr = nf90_put_var(recv_ncid,recv_varid(4),(/sngl(val(1:mesh%nrecv,3))/),(/1,it/),(/mesh%nrecv,1/),(/1,1/))
  !call check2(ierr,'put_var sigma')
  !ierr = nf90_put_var(recv_ncid,recv_varid(5),(/sngl(val(1:mesh%nrecv,4))/),(/1,it/),(/mesh%nrecv,1/),(/1,1/))
  !call check2(ierr,'put_var slip')

  ierr = nf90_sync(recv_ncid)

end subroutine

subroutine recv_io_end(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: ierr
  if (mesh%nrecv==0) return

  ierr = nf90_close(recv_ncid)
  call check2(ierr,'nf90_close recv')

end subroutine

function interp_dist(x,y,z,v,n,x1,y1,z1) result (v1)
  implicit none

  integer :: n,i
  real(kind=rkind),dimension(n) :: x,y,z,v
  real(kind=rkind) :: x1,y1,z1,v1
  real(kind=rkind) :: r1!,r2,r3
  real(kind=rkind) :: w1!,w2,w3,
  real(kind=rkind) :: wsum,vsum

  !r1=sqrt((p(1)-v1(1))**2+(p(2)-v1(2))**2+(p(3)-v1(3))**2)
  !r2=sqrt((p(1)-v2(1))**2+(p(2)-v2(2))**2+(p(3)-v2(3))**2)
  !r3=sqrt((p(1)-v3(1))**2+(p(2)-v3(2))**2+(p(3)-v3(3))**2)
  !
  !w1=1d0/(r1+1d-300)
  !w2=1d0/(r2+1d-300)
  !w3=1d0/(r3+1d-300)
  !
  !cp=(w1*c1+w2*c2+w3*c3)/(w1+w2+w3)

  vsum=0.
  wsum=0.
  do i = 1,n
    r1=sqrt((x(i)-x1)**2+(y(i)-y1)**2+(z(i)-z1)**2)
    w1=1d0/(r1+1d-30)
    vsum=vsum+w1*v(i)
    wsum=wsum+w1
  end do

  v1=vsum/wsum

end function

end module
