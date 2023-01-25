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

module mod_io_grdsurf

  use mod_para,  only : RKIND,                 &
                        Nfp, Np, Nfaces,       &
                        BC_FREE,               &
                        export_grdsurf_displ,  &
                        export_grdsurf_strain, &
                        data_dir
  use mod_types, only : meshvar
  use netcdf

  implicit none

  integer :: grdsurf_ncid
  integer :: grdsurf_varid(20)

  !real(kind=rkind),allocatable,dimension(:,:) :: x,y,z
  real(kind=rkind),allocatable,dimension(:,:) :: PGVx,PGVy,PGVz,PGVh
  real(kind=rkind),allocatable,dimension(:,:) :: Vx,Vy,Vz
  real(kind=rkind),allocatable,dimension(:,:) :: Ux,Uy,Uz
  real(kind=rkind),allocatable,dimension(:,:) :: exx,eyy,ezz,eyz,exz,exy

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

subroutine grdsurf_io_init(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: myrank
  integer :: i, j, ie, is, ierr
  character(len=128) :: filename
  real(kind=rkind),allocatable,dimension(:,:) :: x,y,z!,PGVx,PGVy,PGVz
  !real(kind=rkind),allocatable,dimension(:,:,:) :: vars
  integer :: dimid(20),varid(20)

  if (mesh%nfree_face==0) return

  myrank = mesh%rank

  write(filename,'(a,a,i6.6,a)') trim(data_dir),'/grdsurf_mpi',myrank,'.nc'
  !print*,filename
  allocate(x(Nfp,mesh%nfree_face))
  allocate(y(Nfp,mesh%nfree_face))
  allocate(z(Nfp,mesh%nfree_face))

  allocate(PGVx(Nfp,mesh%nfree_face))
  allocate(PGVy(Nfp,mesh%nfree_face))
  allocate(PGVz(Nfp,mesh%nfree_face))
  allocate(PGVh(Nfp,mesh%nfree_face))
  allocate(Vx(Nfp,mesh%nfree_face))
  allocate(Vy(Nfp,mesh%nfree_face))
  allocate(Vz(Nfp,mesh%nfree_face))
  allocate(Ux(Nfp,mesh%nfree_face))
  allocate(Uy(Nfp,mesh%nfree_face))
  allocate(Uz(Nfp,mesh%nfree_face))
  allocate(exx(Nfp,mesh%nfree_face))
  allocate(eyy(Nfp,mesh%nfree_face))
  allocate(ezz(Nfp,mesh%nfree_face))
  allocate(eyz(Nfp,mesh%nfree_face))
  allocate(exz(Nfp,mesh%nfree_face))
  allocate(exy(Nfp,mesh%nfree_face))
  !allocate(vars(Nfp,mesh%nfree_face,3))
  i = 0
  do ie = 1,mesh%nelem
    do is = 1,Nfaces
      if ( mesh%bctype(is,ie) == BC_FREE) then
        i = i + 1
        x (1:Nfp,i) = mesh%vx(mesh%vmapM(1:Nfp,is,ie))
        y (1:Nfp,i) = mesh%vy(mesh%vmapM(1:Nfp,is,ie))
        z (1:Nfp,i) = mesh%vz(mesh%vmapM(1:Nfp,is,ie))
      endif
    enddo
  enddo
  ! PGV
  do i = 1,Nfp
    do j = 1,mesh%nfree_face
      PGVx(i,j) = -1e30
      PGVy(i,j) = -1e30
      PGVz(i,j) = -1e30
      PGVh(i,j) = -1e30
    enddo
  enddo

  ! start to create
  ierr = nf90_create(filename, NF90_CLOBBER, grdsurf_ncid)
  call check2(ierr,'nf90_create wave')

  ! define dimensions
  ierr = nf90_def_dim(grdsurf_ncid, "Nfp", Nfp, dimid(1))
  call check2(ierr,'def_dim Nfp')
  ierr = nf90_def_dim(grdsurf_ncid, "Nsurf", mesh%nfree_face, dimid(2))
  call check2(ierr,'def_dim Nsurf')
  ierr = nf90_def_dim(grdsurf_ncid, "Nt", nf90_unlimited, dimid(3))
  call check2(ierr,'def_dim Nt')

  ! define variables
  ierr = nf90_def_var(grdsurf_ncid, "x" , NF90_DOUBLE, dimid(1:2), varid(1))
  call check2(ierr,'def_var x')
  ierr = nf90_def_var(grdsurf_ncid, "y" , NF90_DOUBLE, dimid(1:2), varid(2))
  call check2(ierr,'def_var y')
  ierr = nf90_def_var(grdsurf_ncid, "z" , NF90_DOUBLE, dimid(1:2), varid(3))
  call check2(ierr,'def_var z')

  ierr = nf90_def_var(grdsurf_ncid, "time", NF90_DOUBLE, dimid(3), grdsurf_varid(4))
  call check2(ierr,'def_var time')

  ierr = nf90_def_var(grdsurf_ncid, "PGVx" , NF90_DOUBLE, dimid(1:2), grdsurf_varid(14))
  call check2(ierr,'def_var PGVx')
  ierr = nf90_def_var(grdsurf_ncid, "PGVy" , NF90_DOUBLE, dimid(1:2), grdsurf_varid(15))
  call check2(ierr,'def_var PGVy')
  ierr = nf90_def_var(grdsurf_ncid, "PGVz" , NF90_DOUBLE, dimid(1:2), grdsurf_varid(16))
  call check2(ierr,'def_var PGVz')
  ierr = nf90_def_var(grdsurf_ncid, "PGVh" , NF90_DOUBLE, dimid(1:2), grdsurf_varid(17))
  call check2(ierr,'def_var PGVh')

  ierr = nf90_def_var(grdsurf_ncid, "Vx", NF90_FLOAT, dimid(1:3), grdsurf_varid(1))
  call check2(ierr,'def_var Vx')
  ierr = nf90_def_var(grdsurf_ncid, "Vy", NF90_FLOAT, dimid(1:3), grdsurf_varid(2))
  call check2(ierr,'def_var Vy')
  ierr = nf90_def_var(grdsurf_ncid, "Vz", NF90_FLOAT, dimid(1:3), grdsurf_varid(3))
  call check2(ierr,'def_var Vz')

  if (export_grdsurf_displ == 1) then
  ierr = nf90_def_var(grdsurf_ncid, "Ux", NF90_FLOAT, dimid(1:3), grdsurf_varid(5))
  call check2(ierr,'def_var Ux')
  ierr = nf90_def_var(grdsurf_ncid, "Uy", NF90_FLOAT, dimid(1:3), grdsurf_varid(6))
  call check2(ierr,'def_var Uy')
  ierr = nf90_def_var(grdsurf_ncid, "Uz", NF90_FLOAT, dimid(1:3), grdsurf_varid(7))
  call check2(ierr,'def_var Uz')
  end if

  if (export_grdsurf_strain == 1) then
  ierr = nf90_def_var(grdsurf_ncid, "exx", NF90_FLOAT, dimid(1:3), grdsurf_varid(8))
  call check2(ierr,'def_var exx')
  ierr = nf90_def_var(grdsurf_ncid, "eyy", NF90_FLOAT, dimid(1:3), grdsurf_varid(9))
  call check2(ierr,'def_var eyy')
  ierr = nf90_def_var(grdsurf_ncid, "ezz", NF90_FLOAT, dimid(1:3), grdsurf_varid(10))
  call check2(ierr,'def_var ezz')
  ierr = nf90_def_var(grdsurf_ncid, "eyz", NF90_FLOAT, dimid(1:3), grdsurf_varid(11))
  call check2(ierr,'def_var eyz')
  ierr = nf90_def_var(grdsurf_ncid, "exz", NF90_FLOAT, dimid(1:3), grdsurf_varid(12))
  call check2(ierr,'def_var exz')
  ierr = nf90_def_var(grdsurf_ncid, "exy", NF90_FLOAT, dimid(1:3), grdsurf_varid(13))
  call check2(ierr,'def_var exy')
  end if

  ! end of define
  ierr = nf90_enddef(grdsurf_ncid)
  call check2(ierr,'enddef')

  ! put variables
  ierr = nf90_put_var(grdsurf_ncid, varid(1), x)
  call check2(ierr,'put_var x')
  ierr = nf90_put_var(grdsurf_ncid, varid(2), y)
  call check2(ierr,'put_var y')
  ierr = nf90_put_var(grdsurf_ncid, varid(3), z)
  call check2(ierr,'put_var z')

  !!!ierr = nf90_close(grdsurf_ncid)
  deallocate(x,y,z)
  !deallocate(vars)
end subroutine

subroutine grdsurf_io_save(mesh,u,displ,it)
  implicit none
  type(meshvar) :: mesh
  real(kind=rkind) :: u(:,:)
  real(kind=rkind) :: displ(:,:)
  integer :: i, ie, is, it, ierr
  !real,allocatable,dimension(:,:) :: Vx,Vy,Vz
  !real,allocatable,dimension(:,:) :: Ux,Uy,Uz
  !real,allocatable,dimension(:,:) :: exx,eyy,ezz,eyz,exz,exy
  integer,dimension(3) :: start,cnt,stride
  if (mesh%nfree_face==0) return
  i = 0
  do ie = 1,mesh%nelem
    do is = 1,Nfaces
      if ( mesh%bctype(is,ie) == BC_FREE) then
        i = i + 1
        Vx(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),1)/mesh%rho(ie) )
        Vy(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),2)/mesh%rho(ie) )
        Vz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),3)/mesh%rho(ie) )
        exx(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),4))
        eyy(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),5))
        ezz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),6))
        eyz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),4))
        exz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),5))
        exy(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),6))
        Ux(1:Nfp,i) = sngl(displ(mesh%vmapM(1:Nfp,is,ie),1)/mesh%rho(ie) )
        Uy(1:Nfp,i) = sngl(displ(mesh%vmapM(1:Nfp,is,ie),2)/mesh%rho(ie) )
        Uz(1:Nfp,i) = sngl(displ(mesh%vmapM(1:Nfp,is,ie),3)/mesh%rho(ie) )
      endif
    enddo
  enddo
  ! PGV
  !do i = 1,Nfp
  !  do j = 1,mesh%nfree_face
  !    PGVx(i,j) = max(PGVx(i,j),abs(Vx(i,j)))
  !    PGVy(i,j) = max(PGVy(i,j),abs(Vy(i,j)))
  !    PGVz(i,j) = max(PGVz(i,j),abs(Vz(i,j)))
  !  enddo
  !end do

  ierr = nf90_put_var(grdsurf_ncid, grdsurf_varid(14), PGVx)
  call check2(ierr,'put_var PGVx')
  ierr = nf90_put_var(grdsurf_ncid, grdsurf_varid(15), PGVy)
  call check2(ierr,'put_var PGVy')
  ierr = nf90_put_var(grdsurf_ncid, grdsurf_varid(16), PGVz)
  call check2(ierr,'put_var PGVz')
  ierr = nf90_put_var(grdsurf_ncid, grdsurf_varid(17), PGVh)
  call check2(ierr,'put_var PGVh')

  start=(/1,1,it+0/); cnt=(/Nfp,mesh%nfree_face,1/); stride=(/1,1,1/)
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(1),Vx,start,cnt,stride)
  call check2(ierr,'put_var Vx')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(2),Vy,start,cnt,stride)
  call check2(ierr,'put_var Vy')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(3),Vz,start,cnt,stride)
  call check2(ierr,'put_var Vz')

  if (export_grdsurf_displ == 1) then
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(5),Ux,start,cnt,stride)
  call check2(ierr,'put_var Ux')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(6),Uy,start,cnt,stride)
  call check2(ierr,'put_var Uy')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(7),Uz,start,cnt,stride)
  call check2(ierr,'put_var Uz')
  end if

  if (export_grdsurf_strain == 1) then
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(8),exx,start,cnt,stride)
  call check2(ierr,'put_var exx')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(9),eyy,start,cnt,stride)
  call check2(ierr,'put_var eyy')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(10),ezz,start,cnt,stride)
  call check2(ierr,'put_var ezz')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(11),eyz,start,cnt,stride)
  call check2(ierr,'put_var eyz')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(12),exz,start,cnt,stride)
  call check2(ierr,'put_var exz')
  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(13),exy,start,cnt,stride)
  call check2(ierr,'put_var exy')
  end if

  ierr = nf90_put_var(grdsurf_ncid,grdsurf_varid(4), &
      (/mesh%current_time/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var time')

  ierr = nf90_sync(grdsurf_ncid)

  !deallocate(Vx,Vy,Vz)
end subroutine

subroutine update_pgv(mesh,u)
  implicit none
  type(meshvar) :: mesh
  real(kind=rkind) :: u(:,:)
  integer :: i, j, ie, is
  if (mesh%nfree_face==0) return
  i = 0
  do ie = 1,mesh%nelem
    do is = 1,Nfaces
      if ( mesh%bctype(is,ie) == BC_FREE) then
        i = i + 1
        Vx(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),1)/mesh%rho(ie) )
        Vy(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),2)/mesh%rho(ie) )
        Vz(1:Nfp,i) = sngl(u(mesh%vmapM(1:Nfp,is,ie),3)/mesh%rho(ie) )
        do j = 1,Nfp
          PGVx(j,i) = max(PGVx(j,i),abs(Vx(j,i)))
          PGVy(j,i) = max(PGVy(j,i),abs(Vy(j,i)))
          PGVz(j,i) = max(PGVz(j,i),abs(Vz(j,i)))
          PGVh(j,i) = max(PGVh(j,i),sqrt(Vx(j,i)**2+Vy(j,i)**2))
        end do
      endif
    enddo
  enddo
end subroutine

subroutine grdsurf_io_end(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: ierr
  if (mesh%nfree_face==0) return
  ierr = nf90_close(grdsurf_ncid)
  call check2(ierr,'nf90_close wave') 

  deallocate(PGVx,PGVy,PGVz,PGVh)
  deallocate(Vx,Vy,Vz)
  deallocate(Ux,Uy,Uz)
  deallocate(exx,eyy,ezz,eyz,exz,exy)
end subroutine

end module
