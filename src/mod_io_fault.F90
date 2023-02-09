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

module mod_io_fault

  use mod_para,  only : RKIND,           &
                        Nfp, Np, Nfaces, &
                        BC_FAULT,        &
                        data_dir
  use mod_types, only : meshvar
  use netcdf

  implicit none

  integer :: fault_ncid
  integer :: fault_varid(50)

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

subroutine fault_io_init(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: myrank
  integer :: i, ie, ief, is, ierr
  character(len=128) :: filename
  real(kind=rkind),allocatable,dimension(:,:) :: x,y,z,nx,ny,nz
  real(kind=rkind),allocatable,dimension(:,:) :: mx,my,mz,lx,ly,lz
  integer,allocatable,dimension(:) :: fault_id
  real(kind=rkind),allocatable,dimension(:,:,:) :: fltvars
  integer :: dimid(50),varid(50)

  if (mesh%nfault_face==0) return

  myrank = mesh%rank

  write(filename,'(a,a,i6.6,a)') trim(data_dir),'/fault_mpi',myrank,'.nc'
  !print*,filename
  allocate(fault_id(mesh%nfault_face))
  allocate(x(Nfp,mesh%nfault_face))
  allocate(y(Nfp,mesh%nfault_face))
  allocate(z(Nfp,mesh%nfault_face))
  allocate(nx(Nfp,mesh%nfault_face))
  allocate(ny(Nfp,mesh%nfault_face))
  allocate(nz(Nfp,mesh%nfault_face))
  allocate(mx(Nfp,mesh%nfault_face))
  allocate(my(Nfp,mesh%nfault_face))
  allocate(mz(Nfp,mesh%nfault_face))
  allocate(lx(Nfp,mesh%nfault_face))
  allocate(ly(Nfp,mesh%nfault_face))
  allocate(lz(Nfp,mesh%nfault_face))
  allocate(fltvars(Nfp,mesh%nfault_face,7))
  i = 0
  !do ie = 1,mesh%nelem
  do ief = 1,mesh%nfault_elem
    ie = mesh%fault2wave(ief)
    do is = 1,Nfaces
      if ( mesh%bctype(is,ie) >= BC_FAULT) then
        i = i + 1
        x (1:Nfp,i) = mesh%vx(mesh%vmapM(1:Nfp,is,ie))
        y (1:Nfp,i) = mesh%vy(mesh%vmapM(1:Nfp,is,ie))
        z (1:Nfp,i) = mesh%vz(mesh%vmapM(1:Nfp,is,ie))
        nx(1:Nfp,i) = mesh%nx(((is-1)*Nfp+1):(is*Nfp),ie)
        ny(1:Nfp,i) = mesh%ny(((is-1)*Nfp+1):(is*Nfp),ie)
        nz(1:Nfp,i) = mesh%nz(((is-1)*Nfp+1):(is*Nfp),ie)
        mx(1:Nfp,i) = mesh%mx(((is-1)*Nfp+1):(is*Nfp),ie)
        my(1:Nfp,i) = mesh%my(((is-1)*Nfp+1):(is*Nfp),ie)
        mz(1:Nfp,i) = mesh%mz(((is-1)*Nfp+1):(is*Nfp),ie)
        lx(1:Nfp,i) = mesh%lx(((is-1)*Nfp+1):(is*Nfp),ie)
        ly(1:Nfp,i) = mesh%ly(((is-1)*Nfp+1):(is*Nfp),ie)
        lz(1:Nfp,i) = mesh%lz(((is-1)*Nfp+1):(is*Nfp),ie)
        fault_id(i) = mesh%bctype(is,ie)-BC_FAULT+1
        fltvars(1:Nfp,i,1) = mesh%mu_s (1:Nfp,is,ief)
        fltvars(1:Nfp,i,2) = mesh%mu_d (1:Nfp,is,ief)
        fltvars(1:Nfp,i,3) = mesh%Dc   (1:Nfp,is,ief)
        fltvars(1:Nfp,i,4) = mesh%C0   (1:Nfp,is,ief)
        fltvars(1:Nfp,i,5) = mesh%tau0n(1:Nfp,is,ief)
        fltvars(1:Nfp,i,6) = mesh%tau0m(1:Nfp,is,ief)
        fltvars(1:Nfp,i,7) = mesh%tau0l(1:Nfp,is,ief)
      endif
    enddo
  enddo

  ! start to create
  ierr = nf90_create(filename, NF90_CLOBBER, fault_ncid)
  call check2(ierr,'nf90_create fault')

  ! define dimensions
  ierr = nf90_def_dim(fault_ncid, "Nfp", Nfp, dimid(1))
  call check2(ierr,'def_dim Nfp')
  ierr = nf90_def_dim(fault_ncid, "Nflt", mesh%nfault_face, dimid(2))
  call check2(ierr,'def_dim Nflt')
  ierr = nf90_def_dim(fault_ncid, "Nt", nf90_unlimited, dimid(3))
  call check2(ierr,'def_dim Nt')

  ! define variables
  ierr = nf90_def_var(fault_ncid, "x" , NF90_DOUBLE, dimid(1:2), varid(1))
  call check2(ierr,'def_var x')
  ierr = nf90_def_var(fault_ncid, "y" , NF90_DOUBLE, dimid(1:2), varid(2))
  call check2(ierr,'def_var y')
  ierr = nf90_def_var(fault_ncid, "z" , NF90_DOUBLE, dimid(1:2), varid(3))
  call check2(ierr,'def_var z')
  ierr = nf90_def_var(fault_ncid, "nx", NF90_DOUBLE, dimid(1:2), varid(4))
  call check2(ierr,'def_var nx')
  ierr = nf90_def_var(fault_ncid, "ny", NF90_DOUBLE, dimid(1:2), varid(5))
  call check2(ierr,'def_var ny')
  ierr = nf90_def_var(fault_ncid, "nz", NF90_DOUBLE, dimid(1:2), varid(6))
  call check2(ierr,'def_var nz')
  ierr = nf90_def_var(fault_ncid, "mx", NF90_DOUBLE, dimid(1:2), varid(15))
  call check2(ierr,'def_var mx')
  ierr = nf90_def_var(fault_ncid, "my", NF90_DOUBLE, dimid(1:2), varid(16))
  call check2(ierr,'def_var my')
  ierr = nf90_def_var(fault_ncid, "mz", NF90_DOUBLE, dimid(1:2), varid(17))
  call check2(ierr,'def_var mz')
  ierr = nf90_def_var(fault_ncid, "lx", NF90_DOUBLE, dimid(1:2), varid(18))
  call check2(ierr,'def_var lx')
  ierr = nf90_def_var(fault_ncid, "ly", NF90_DOUBLE, dimid(1:2), varid(19))
  call check2(ierr,'def_var ly')
  ierr = nf90_def_var(fault_ncid, "lz", NF90_DOUBLE, dimid(1:2), varid(20))
  call check2(ierr,'def_var lz')
  ierr = nf90_def_var(fault_ncid, "mu_s", NF90_DOUBLE, dimid(1:2), varid(7))
  call check2(ierr,'def_var mu_s')
  ierr = nf90_def_var(fault_ncid, "mu_d", NF90_DOUBLE, dimid(1:2), varid(8))
  call check2(ierr,'def_var mu_d')
  ierr = nf90_def_var(fault_ncid, "Dc"  , NF90_DOUBLE, dimid(1:2), varid(9))
  call check2(ierr,'def_var Dc')
  ierr = nf90_def_var(fault_ncid, "C0"  , NF90_DOUBLE, dimid(1:2), varid(10))
  call check2(ierr,'def_var C0')
  ierr = nf90_def_var(fault_ncid, "tau0n", NF90_DOUBLE, dimid(1:2), varid(11))
  call check2(ierr,'def_var tau0n')
  ierr = nf90_def_var(fault_ncid, "tau0m", NF90_DOUBLE, dimid(1:2), varid(12))
  call check2(ierr,'def_var tau0m')
  ierr = nf90_def_var(fault_ncid, "tau0l", NF90_DOUBLE, dimid(1:2), varid(13))
  call check2(ierr,'def_var tau0l')
  ierr = nf90_def_var(fault_ncid, "fault_id", NF90_INT, dimid(2), varid(14))
  call check2(ierr,'def_var fault_id')
  ierr = nf90_def_var(fault_ncid, "ruptime", NF90_DOUBLE, dimid(1:2), fault_varid(4))
  call check2(ierr,'def_var ruptime')
  ierr = nf90_def_var(fault_ncid, "peakrate", NF90_DOUBLE, dimid(1:2), fault_varid(16))
  call check2(ierr,'def_var peakrate')

  ierr = nf90_def_var(fault_ncid, "time" , NF90_DOUBLE, dimid(3), fault_varid(8))
  call check2(ierr,'def_var time')

  ierr = nf90_def_var(fault_ncid, "rate",   NF90_FLOAT, dimid(1:3), fault_varid(1))
  call check2(ierr,'def_var rate')
  ierr = nf90_def_var(fault_ncid, "ratem",  NF90_FLOAT, dimid(1:3), fault_varid(5))
  call check2(ierr,'def_var ratem')
  ierr = nf90_def_var(fault_ncid, "ratel",  NF90_FLOAT, dimid(1:3), fault_varid(6))
  call check2(ierr,'def_var ratel')
  ierr = nf90_def_var(fault_ncid, "slip",  NF90_FLOAT, dimid(1:3), fault_varid(7))
  call check2(ierr,'def_var slip')
  ierr = nf90_def_var(fault_ncid, "slipm",  NF90_FLOAT, dimid(1:3), fault_varid(14))
  call check2(ierr,'def_var slipm')
  ierr = nf90_def_var(fault_ncid, "slipl",  NF90_FLOAT, dimid(1:3), fault_varid(15))
  call check2(ierr,'def_var slipl')
  ierr = nf90_def_var(fault_ncid, "tau", NF90_FLOAT, dimid(1:3), fault_varid(2))
  call check2(ierr,'def_var tau')
  ierr = nf90_def_var(fault_ncid, "taum",  NF90_FLOAT, dimid(1:3), fault_varid(12))
  call check2(ierr,'def_var taum')
  ierr = nf90_def_var(fault_ncid, "taul",  NF90_FLOAT, dimid(1:3), fault_varid(13))
  call check2(ierr,'def_var taul')
  ierr = nf90_def_var(fault_ncid, "sigma",  NF90_FLOAT, dimid(1:3), fault_varid(3))
  call check2(ierr,'def_var sigma')
  ! rate state
  ierr = nf90_def_var(fault_ncid, "state",  NF90_FLOAT, dimid(1:3), fault_varid(9))
  call check2(ierr,'def_var state')
  ! thermal pressurization
  ierr = nf90_def_var(fault_ncid, "TP_T",  NF90_DOUBLE, dimid(1:3), fault_varid(10))
  call check2(ierr,'def_var TP_T')
  ierr = nf90_def_var(fault_ncid, "TP_P",  NF90_DOUBLE, dimid(1:3), fault_varid(11))
  call check2(ierr,'def_var TP_P')

  ! end of define
  ierr = nf90_enddef(fault_ncid)
  call check2(ierr,'enddef')

  ! put variables
  ierr = nf90_put_var(fault_ncid, varid(1), x)
  call check2(ierr,'put_var x')
  ierr = nf90_put_var(fault_ncid, varid(2), y)
  call check2(ierr,'put_var y')
  ierr = nf90_put_var(fault_ncid, varid(3), z)
  call check2(ierr,'put_var z')
  ierr = nf90_put_var(fault_ncid, varid(4), nx)
  call check2(ierr,'put_var nx')
  ierr = nf90_put_var(fault_ncid, varid(5), ny)
  call check2(ierr,'put_var ny')
  ierr = nf90_put_var(fault_ncid, varid(6), nz)
  call check2(ierr,'put_var nz')
  ierr = nf90_put_var(fault_ncid, varid(15), mx)
  call check2(ierr,'put_var mx')
  ierr = nf90_put_var(fault_ncid, varid(16), my)
  call check2(ierr,'put_var my')
  ierr = nf90_put_var(fault_ncid, varid(17), mz)
  call check2(ierr,'put_var mz')
  ierr = nf90_put_var(fault_ncid, varid(18), lx)
  call check2(ierr,'put_var lx')
  ierr = nf90_put_var(fault_ncid, varid(19), ly)
  call check2(ierr,'put_var ly')
  ierr = nf90_put_var(fault_ncid, varid(20), lz)
  call check2(ierr,'put_var lz')
  ierr = nf90_put_var(fault_ncid, varid(7), fltvars(:,:,1))
  call check2(ierr,'put_var mu_s')
  ierr = nf90_put_var(fault_ncid, varid(8), fltvars(:,:,2))
  call check2(ierr,'put_var mu_d')
  ierr = nf90_put_var(fault_ncid, varid(9), fltvars(:,:,3))
  call check2(ierr,'put_var Dc')
  ierr = nf90_put_var(fault_ncid, varid(10), fltvars(:,:,4))
  call check2(ierr,'put_var C0')
  ierr = nf90_put_var(fault_ncid, varid(11), fltvars(:,:,5))
  call check2(ierr,'put_var tau0n')
  ierr = nf90_put_var(fault_ncid, varid(12), fltvars(:,:,6))
  call check2(ierr,'put_var tau0m')
  ierr = nf90_put_var(fault_ncid, varid(13), fltvars(:,:,7))
  call check2(ierr,'put_var tau0l')
  ierr = nf90_put_var(fault_ncid, varid(14), fault_id)
  call check2(ierr,'put_var fault_id')

  !!!ierr = nf90_close(fault_ncid)
  deallocate(x,y,z,nx,ny,nz)
  deallocate(fltvars)
  deallocate(fault_id)
end subroutine

subroutine fault_io_save(mesh,it)
  implicit none
  type(meshvar) :: mesh
  integer :: i, ie, ief, is, it, ierr
  real,allocatable,dimension(:,:) :: rate,stress,sigma,ruptime,rate1,rate2,slip,peakrate
  real,allocatable,dimension(:,:) :: stress1,stress2
  real,allocatable,dimension(:,:) :: slip1,slip2
  ! rate state
  real,allocatable,dimension(:,:) :: state
  ! thermal pressurization
  real*8,allocatable,dimension(:,:) :: TP_T,TP_P
  integer,dimension(3) :: start,cnt,stride
  if (mesh%nfault_face==0) return
  allocate(rate   (Nfp,mesh%nfault_face))
  allocate(rate1  (Nfp,mesh%nfault_face))
  allocate(rate2  (Nfp,mesh%nfault_face))
  allocate(slip   (Nfp,mesh%nfault_face))
  allocate(slip1  (Nfp,mesh%nfault_face))
  allocate(slip2  (Nfp,mesh%nfault_face))
  allocate(stress (Nfp,mesh%nfault_face))
  allocate(stress1(Nfp,mesh%nfault_face))
  allocate(stress2(Nfp,mesh%nfault_face))
  allocate(sigma  (Nfp,mesh%nfault_face))
  allocate(ruptime(Nfp,mesh%nfault_face))
  allocate(peakrate(Nfp,mesh%nfault_face))
  ! rate state
  allocate(state(Nfp,mesh%nfault_face))
  ! thermal pressurization
  allocate(TP_T(Nfp,mesh%nfault_face))
  allocate(TP_P(Nfp,mesh%nfault_face))

  i = 0
  !do ie = 1,mesh%nelem
  do ief = 1,mesh%nfault_elem
    ie = mesh%fault2wave(ief)
    do is = 1,Nfaces
      if ( mesh%bctype(is,ie) >= BC_FAULT) then
        i = i + 1
        rate   (1:Nfp,i) = sngl(mesh%sliprate (1:Nfp,is,ief))
        stress (1:Nfp,i) = sngl(mesh%stress   (1:Nfp,is,ief))
        stress1(1:Nfp,i) = sngl(mesh%stress1  (1:Nfp,is,ief))
        stress2(1:Nfp,i) = sngl(mesh%stress2  (1:Nfp,is,ief))
        sigma  (1:Nfp,i) = sngl(mesh%sigma    (1:Nfp,is,ief))
        rate1  (1:Nfp,i) = sngl(mesh%sliprate1(1:Nfp,is,ief))
        rate2  (1:Nfp,i) = sngl(mesh%sliprate2(1:Nfp,is,ief))
        slip   (1:Nfp,i) = sngl(mesh%slip     (1:Nfp,is,ief))
        slip1  (1:Nfp,i) = sngl(mesh%slip1    (1:Nfp,is,ief))
        slip2  (1:Nfp,i) = sngl(mesh%slip2    (1:Nfp,is,ief))
        ruptime(1:Nfp,i) = sngl(mesh%ruptime  (1:Nfp,is,ief))
        peakrate(1:Nfp,i) = sngl(mesh%peakrate(1:Nfp,is,ief))
        ! rate state
        state(1:Nfp,i) = sngl(mesh%state(1:Nfp,is,ief))
        ! thermal pressurization
        TP_T(1:Nfp,i) = (mesh%TP_T(1:Nfp,is,ief))
        TP_P(1:Nfp,i) = (mesh%TP_P(1:Nfp,is,ief))

      endif
    enddo
  enddo
  start=(/1,1,it+0/); cnt=(/Nfp,mesh%nfault_face,1/); stride=(/1,1,1/)
  ierr = nf90_put_var(fault_ncid,fault_varid(1),rate,start,cnt,stride)
  call check2(ierr,'put_var rate')
  ierr = nf90_put_var(fault_ncid,fault_varid(2),stress,start,cnt,stride)
  call check2(ierr,'put_var tau')
  ierr = nf90_put_var(fault_ncid,fault_varid(3),sigma,start,cnt,stride)
  call check2(ierr,'put_var sigma')
  ierr = nf90_put_var(fault_ncid,fault_varid(4),ruptime)
  call check2(ierr,'put_var ruptime')
  ierr = nf90_put_var(fault_ncid,fault_varid(16),peakrate)
  call check2(ierr,'put_var peakrate')
  ierr = nf90_put_var(fault_ncid,fault_varid(5),rate1,start,cnt,stride)
  call check2(ierr,'put_var ratem')
  ierr = nf90_put_var(fault_ncid,fault_varid(6),rate2,start,cnt,stride)
  call check2(ierr,'put_var raten')
  ierr = nf90_put_var(fault_ncid,fault_varid(12),stress1,start,cnt,stride)
  call check2(ierr,'put_var taum')
  ierr = nf90_put_var(fault_ncid,fault_varid(13),stress2,start,cnt,stride)
  call check2(ierr,'put_var taul')
  ierr = nf90_put_var(fault_ncid,fault_varid(7),slip,start,cnt,stride)
  call check2(ierr,'put_var slip')
  ierr = nf90_put_var(fault_ncid,fault_varid(14),slip1,start,cnt,stride)
  call check2(ierr,'put_var slip1')
  ierr = nf90_put_var(fault_ncid,fault_varid(15),slip2,start,cnt,stride)
  call check2(ierr,'put_var slip2')
  ! rate state
  ierr = nf90_put_var(fault_ncid,fault_varid(9),state,start,cnt,stride)
  call check2(ierr,'put_var state')
  ! thermal pressurization
  ierr = nf90_put_var(fault_ncid,fault_varid(10),TP_T,start,cnt,stride)
  call check2(ierr,'put_var TP_T')
  ierr = nf90_put_var(fault_ncid,fault_varid(11),TP_P,start,cnt,stride)
  call check2(ierr,'put_var TP_P')

  ierr = nf90_put_var(fault_ncid,fault_varid(8), &
      (/mesh%current_time/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var time')

  ierr = nf90_sync(fault_ncid)

  deallocate(rate,stress,sigma,ruptime,slip)
  deallocate(rate1,rate2)
  deallocate(state)
  deallocate(TP_T,TP_P)
  deallocate(peakrate)
end subroutine

subroutine fault_io_end(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: ierr
  if (mesh%nfault_face==0) return
  ierr = nf90_close(fault_ncid)
  call check2(ierr,'nf90_close fault')
end subroutine

end module
