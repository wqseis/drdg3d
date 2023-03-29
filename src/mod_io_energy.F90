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

module mod_io_energy

  use mod_para,  only : RKIND,           &
                        Nfp, Np, Nfaces, &
                        BC_FAULT,        &
                        data_dir
  use mod_types, only : meshvar
  use mod_energy, only : cal_energy, cal_moment
  use netcdf

  implicit none

  integer :: en_ncid
  integer :: en_varid(5)

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

subroutine energy_io_init(mesh)
  implicit none
  type(meshvar) :: mesh
  character(len=128) :: filename
  integer :: dimid(2)
  integer :: ierr

  if (mesh%rank .ne. 0) return

  write(filename,*) trim(data_dir),'/energy.nc'

  ! start to create
  ierr = nf90_create(filename, NF90_CLOBBER, en_ncid)
  call check2(ierr,'nf90_create energy')

  ! define dimensions
  ierr = nf90_def_dim(en_ncid, "Nt", nf90_unlimited, dimid(1))
  call check2(ierr,'def_dim Nt')

  ! define variables

  ierr = nf90_def_var(en_ncid, "time", NF90_DOUBLE, dimid(1), en_varid(1))
  call check2(ierr,'def_var time')
  ierr = nf90_def_var(en_ncid, "KinematicEnergy", NF90_DOUBLE, dimid(1), en_varid(2))
  call check2(ierr,'def_var Ev')
  ierr = nf90_def_var(en_ncid, "StrainEnergy", NF90_DOUBLE, dimid(1), en_varid(3))
  call check2(ierr,'def_var Es')
  ierr = nf90_def_var(en_ncid, "Moment", NF90_DOUBLE, dimid(1), en_varid(4))
  call check2(ierr,'def_var Em')
  ierr = nf90_def_var(en_ncid, "MomentRate", NF90_DOUBLE, dimid(1), en_varid(5))
  call check2(ierr,'def_var Er')

  ! end of define
  ierr = nf90_enddef(en_ncid)
  call check2(ierr,'enddef')

end subroutine

subroutine energy_io_save(mesh,u,it)
  implicit none
  type(meshvar) :: mesh
  real(kind=RKIND),dimension(:,:) :: u
  real(kind=rkind) :: Ev,Es,Mmt,Mmtr
  integer :: ierr, it

  call cal_energy(mesh,u,Ev,Es)
  call cal_moment(mesh,Mmt,Mmtr)

  if (mesh%rank == 0) then

  ierr = nf90_put_var(en_ncid,en_varid(1),(/mesh%current_time/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var time')
  ierr = nf90_put_var(en_ncid,en_varid(2),(/Ev/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var Ev')
  ierr = nf90_put_var(en_ncid,en_varid(3),(/Es/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var Es')
  ierr = nf90_put_var(en_ncid,en_varid(4),(/Mmt/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var Em')
  ierr = nf90_put_var(en_ncid,en_varid(5),(/Mmtr/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var Er')

  ierr = nf90_sync(en_ncid)

  end if

end subroutine

subroutine energy_io_end(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: ierr
  if (mesh%rank .ne. 0) return
  ierr = nf90_close(en_ncid)
  call check2(ierr,'nf90_close energy')
end subroutine

end module
