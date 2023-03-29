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

module mod_mesh

  use netcdf
  use mod_para,  only : RKIND,                    &
                        Nfaces,                   &
                        mesh_dir,                 &
                        export_media
  use mod_types, only : meshvar
  use mod_vtk,   only : writeVtkTetraMeshRealdata

  implicit none

  ! public :: meshvar
  ! 
  ! type :: meshvar
  !   integer :: nelem,ncoord
  !   integer :: nglob
  ! 
  !   integer :: mpi_nn
  !   integer :: mpi_ne
  !   integer :: mpi_nemax
  !   integer :: pinterfaces
  !   integer :: nrecv
  !   integer :: body_nrecv
  !   integer :: nfault_elem
  !   integer :: irk,nrk
  ! 
  !   real(kind=RKIND) :: dtfactor,deltat
  !   real(kind=RKIND) :: current_time
  ! 
  !   real(kind=RKIND),dimension(:,:),pointer :: coord=>null()
  !   integer,dimension(:,:),pointer :: elem => null()
  !   integer,dimension(:,:),pointer :: neigh => null()
  !   integer,dimension(:,:),pointer :: face => null()
  !   integer,dimension(:,:),pointer :: direction => null()
  !   integer,dimension(:,:),pointer :: bctype => null()
  !   integer,dimension(:),pointer :: elemtype => null()
  !   !real(kind=RKIND),dimension(:,:),allocatable :: coord
  !   !integer,dimension(:,:),allocatable :: elem
  !   !integer,dimension(:,:),allocatable :: neigh
  !   !integer,dimension(:,:),allocatable :: face
  !   !integer,dimension(:,:),allocatable :: direction
  !   !integer,dimension(:,:),allocatable :: bctype
  !   !integer,dimension(:),allocatable :: elemtype
  !   !integer :: vmapM(Nfp,Nfaces)
  !   integer,dimension(:,:,:),pointer :: vmapM => null()
  !   integer,dimension(:,:,:),pointer :: vmapP => null()
  !   !integer :: flipped_index(Nfp,6)
  !   integer,dimension(:,:),pointer :: flipped_index => null()
  !   ! mpi
  !   integer,dimension(:),pointer :: mpi_neighbor => null()
  !   integer,dimension(:,:,:),pointer :: mpi_connection => null()
  !   integer,dimension(:,:),pointer :: mpi_ibool => null()
  !   integer,dimension(:,:,:),pointer :: mpi_interface => null()
  !   real(kind=rkind),dimension(:,:),pointer :: mpi_vp => null()
  !   real(kind=rkind),dimension(:,:),pointer :: mpi_vs => null()
  !   real(kind=rkind),dimension(:,:),pointer :: mpi_rho => null()
  !   !integer,dimension(:),allocatable :: mpi_neighbor
  !   !integer,dimension(:,:,:),allocatable :: mpi_connection
  !   !integer,dimension(:,:),allocatable :: mpi_ibool
  !   !integer,dimension(:,:,:),allocatable :: mpi_interface
  ! 
  !   ! fault
  !   real(kind=rkind),dimension(:,:,:),pointer :: tau0n => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: tau0m => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: tau0l => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: stress => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: sigma => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: slip => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: mslip => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: tslip => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: sliprate => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: sliprate1 => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: sliprate2 => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: mu_s => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: mu_d => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: Dc => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: C0 => null()
  !   real(kind=rkind),dimension(:,:,:),pointer :: ruptime => null()
  ! 
  !   ! matrices
  !   !real(kind=rkind),dimension(Np,Np) :: Vdm,invVdm,Mass,Dr,Ds,Dt!,Filter
  !   !integer :: Fmask(Nfp,Nfaces)
  !   !real(kind=rkind) :: LIFT(Np,Nfaces*Nfp)
  !   integer,dimension(:,:),pointer :: Fmask => null()
  !   real(kind=rkind),dimension(:,:),pointer :: Vdm => null()
  !   real(kind=rkind),dimension(:,:),pointer :: invVdm => null()
  !   real(kind=rkind),dimension(:,:),pointer :: Mass => null()
  !   real(kind=rkind),dimension(:,:),pointer :: Dr => null()
  !   real(kind=rkind),dimension(:,:),pointer :: Ds => null()
  !   real(kind=rkind),dimension(:,:),pointer :: Dt => null()
  !   real(kind=rkind),dimension(:,:),pointer :: LIFT => null()
  ! 
  ! 
  !   real(kind=rkind),dimension(:),pointer :: r => null()
  !   real(kind=rkind),dimension(:),pointer :: s => null()
  !   real(kind=rkind),dimension(:),pointer :: t => null()
  !   real(kind=rkind),dimension(:),pointer :: vx => null()
  !   real(kind=rkind),dimension(:),pointer :: vy => null()
  !   real(kind=rkind),dimension(:),pointer :: vz => null()
  !   real(kind=rkind),dimension(:),pointer :: rx => null()
  !   real(kind=rkind),dimension(:),pointer :: ry => null()
  !   real(kind=rkind),dimension(:),pointer :: rz => null()
  !   real(kind=rkind),dimension(:),pointer :: sx => null()
  !   real(kind=rkind),dimension(:),pointer :: sy => null()
  !   real(kind=rkind),dimension(:),pointer :: sz => null()
  !   real(kind=rkind),dimension(:),pointer :: tx => null()
  !   real(kind=rkind),dimension(:),pointer :: ty => null()
  !   real(kind=rkind),dimension(:),pointer :: tz => null()
  !   real(kind=rkind),dimension(:),pointer :: jac => null()
  !   real(kind=rkind),dimension(:,:),pointer :: nx => null()
  !   real(kind=rkind),dimension(:,:),pointer :: ny => null()
  !   real(kind=rkind),dimension(:,:),pointer :: nz => null()
  !   real(kind=rkind),dimension(:,:),pointer :: mx => null()
  !   real(kind=rkind),dimension(:,:),pointer :: my => null()
  !   real(kind=rkind),dimension(:,:),pointer :: mz => null()
  !   real(kind=rkind),dimension(:,:),pointer :: lx => null()
  !   real(kind=rkind),dimension(:,:),pointer :: ly => null()
  !   real(kind=rkind),dimension(:,:),pointer :: lz => null()
  !   real(kind=rkind),dimension(:,:),pointer :: sJ => null()
  !   real(kind=rkind),dimension(:,:),pointer :: Fscale => null()
  ! 
  !   ! media
  !   real(kind=rkind),dimension(:),pointer :: vp => null()
  !   real(kind=rkind),dimension(:),pointer :: vs => null()
  !   real(kind=rkind),dimension(:),pointer :: rho => null()
  !   real(kind=rkind),dimension(:),pointer :: zp => null()
  !   real(kind=rkind),dimension(:),pointer :: zs => null()
  !   real(kind=rkind),dimension(:),pointer :: lam => null()
  !   real(kind=rkind),dimension(:),pointer :: mu => null()
  !   !real(kind=rkind),dimension(:),allocatable :: vp
  !   !real(kind=rkind),dimension(:),allocatable :: vs
  !   !real(kind=rkind),dimension(:),allocatable :: rho
  !   !real(kind=rkind),dimension(:),allocatable :: zp
  !   !real(kind=rkind),dimension(:),allocatable :: zs
  !   !real(kind=rkind),dimension(:),allocatable :: lam
  !   !real(kind=rkind),dimension(:),allocatable :: mu
  ! !@
  ! !@  ! fault recvs
  ! !@  integer,dimension(:),pointer :: recv_fid => null()
  ! !@  integer,dimension(:),pointer :: recv_i => null()
  ! !@  integer,dimension(:),pointer :: recv_ie => null()
  ! !@  real(kind=rkind),dimension(:),pointer :: recv_refx => null()
  ! !@  real(kind=rkind),dimension(:,:,:),pointer :: recv_buffer => null()
  ! !@  real(kind=rkind),dimension(:,:,:),pointer :: fault_buffer => null()
  ! !@
  ! !@  real(kind=rkind),dimension(:,:,:),pointer :: surface_buffer => null()
  ! !@
  ! !@  ! body recvs
  ! !@  integer,dimension(:),pointer :: body_recv_fid => null()
  ! !@  integer,dimension(:),pointer :: body_recv_i => null()
  ! !@  integer,dimension(:),pointer :: body_recv_j => null()
  ! !@  integer,dimension(:),pointer :: body_recv_ie => null()
  ! !@  real(kind=rkind),dimension(:),pointer :: body_recv_refx => null()
  ! !@  real(kind=rkind),dimension(:),pointer :: body_recv_refy => null()
  ! end type

contains

subroutine check2(status,msg)
  implicit none
  integer, intent ( in) :: status
  character(len=*) :: msg

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))!,' in file ',__FILE__,' line ',__LINE__
    print *, trim(msg)
    stop 110
  end if
end subroutine

subroutine readMeshVar(this)
  implicit none
  type(meshvar),intent(inout) :: this
  integer :: myrank,nproc

  !integer :: ios
  character(len=80) :: filename,namestr
  integer :: ierr,ncid,dimid,varid
  !integer :: i,j,ie,is!,k!,is
  !real(kind=rkind),dimension(NGLL) :: xface,yface

  myrank = this%rank
  nproc = this%nproc

  !write(filename,'(a,i6.6)') 'data/meshVar',myrank
  !write(filename,'(a,i6.6,a)') 'data/meshVar',myrank,'.nc'
  write(filename,'(a,a,i6.6,a)') trim(mesh_dir),'/meshVar',myrank,'.nc'
  !print*,trim(filename)
  !open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown', iostat = ios)
!  open(unit=27,file=trim(filename),status='old',iostat=ios)
!  if (ios /= 0) then
!      print*,'could not open: '//trim(filename)
!  end if

!  read(27,*) this%ncoord

!  allocate(this%coord(3,this%ncoord))

  !read(27,*) this%coord

  !read(27,*) this%nelem

  !allocate(this%elem     (     4,this%nelem))
!@  allocate(this%neigh    (Nfaces,this%nelem))
!@  allocate(this%face     (Nfaces,this%nelem))
!@  allocate(this%direction(Nfaces,this%nelem))
!@  allocate(this%bctype   (Nfaces,this%nelem))
!@  allocate(this%elemtype (       this%nelem))
!@  ! media
!@  allocate(this%rho(this%nelem))
!@  allocate(this%vp (this%nelem))
!@  allocate(this%vs (this%nelem))
!@  allocate(this%zp (this%nelem))
!@  allocate(this%zs (this%nelem))
!@  allocate(this%mu (this%nelem))
!@  allocate(this%lam(this%nelem))

  !this%elemtype = ELEM_FLUID
  !this%elemtype = ELEM_SOLID

!  read(27,*) this%elem
!@  read(27,*) this%neigh
!@  read(27,*) this%face
!@  read(27,*) this%direction
!@  read(27,*) this%bctype
!@  read(27,*) this%elemtype
!@  read(27,*) this%rho
!@  read(27,*) this%vp
!@  read(27,*) this%vs

  !do i = 1,this%nelem
  !write(*,*) this%bctype(:,i)
  !end do
!if (.false.) then
!@  read(27,*) this%mpi_nn
!@  read(27,*) this%mpi_ne
!@  read(27,*) this%mpi_nemax
!@
!@  allocate(this%mpi_neighbor(this%mpi_nn))
!@  allocate(this%mpi_connection(this%mpi_nn,this%mpi_ne,2))
!@  allocate(this%mpi_ibool(Nfaces,this%mpi_nemax))
!@  allocate(this%mpi_interface(4,Nfaces,this%nelem))
!@
!@  read(27,*) this%mpi_neighbor
!@  read(27,*) this%mpi_connection
!@  read(27,*) this%mpi_ibool
!@  read(27,*) this%mpi_interface
!@
!@
!@  read(27,*) this%nrecv
!@  if (this%nrecv > 0) then
!@      allocate(this%recv_fid   (this%nrecv))
!@      allocate(this%recv_i     (this%nrecv))
!@      allocate(this%recv_ie    (this%nrecv))
!@      allocate(this%recv_refx  (this%nrecv))
!@      allocate(this%recv_buffer(this%nrecv,NGLL,10))
!@
!@      read(27,*) this%recv_fid
!@      read(27,*) this%recv_i
!@      read(27,*) this%recv_ie
!@      read(27,*) this%recv_refx
!@  end if
!@
!@  read(27,*) this%body_nrecv
!@  ! print*,'rank=',myrank,'body_nrecv=',this%body_nrecv
!@  if (this%body_nrecv > 0) then
!@      allocate(this%body_recv_fid (this%body_nrecv))
!@      allocate(this%body_recv_i   (this%body_nrecv))
!@      allocate(this%body_recv_j   (this%body_nrecv))
!@      allocate(this%body_recv_ie  (this%body_nrecv))
!@      allocate(this%body_recv_refx(this%body_nrecv))
!@      allocate(this%body_recv_refy(this%body_nrecv))
!@
!@      read(27,*) this%body_recv_fid
!@      read(27,*) this%body_recv_i
!@      read(27,*) this%body_recv_j
!@      read(27,*) this%body_recv_ie
!@      read(27,*) this%body_recv_refx
!@      read(27,*) this%body_recv_refy
!@  end if
!@
!end if

!  close(27)

!@  if (.false.) then
!@
!@      print*,'ncoord=',this%ncoord
!@      do j = 1,this%ncoord
!@      print*,(this%coord(i,j),i=1,3)
!@      enddo
!@
!@      print*,'nelem=',this%nelem
!@      do j = 1,this%nelem
!@      print*,(this%elem(i,j),i=1,4)
!@      enddo
!@
!@      print*,'neigh=',this%nelem
!@      do j = 1,this%nelem
!@      print*,(this%neigh(i,j),i=1,Nfaces)
!@      enddo
!@
!@      print*,'face=',this%nelem
!@      do j = 1,this%nelem
!@      print*,(this%face(i,j),i=1,Nfaces)
!@      enddo
!@
!@      stop 2
!@  endif

  !k = 0
  !do i = 1,this%nelem
  !    do j = 1,Nfaces
  !        if (this%bctype(j,i) >= BC_FAULT) k=k+1
  !    end do
  !end do
  !this%nfault_face = k
!@  ! check direction
!@  do ie = 1,this%nelem
!@    do is = 1,Nfaces
!@      if( this%direction(is,ie) < 1 .or. this%direction(is,ie) > 6 ) then
!@        print*,'error at direction, rank=',myrank,'ie=',ie,'is=',is,'dire=',this%direction(is,ie)
!@        stop 111
!@      end if
!@    end do
!@  end do
  ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid)
  call check2(ierr, 'nf90_open '//trim(filename))

  ierr = nf90_inq_dimid(ncid,'Nelem',dimid)
  call check2(ierr, 'inq_dimid Nelem')

  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%nelem)
  call check2(ierr, 'inquire_dimension Nelem')

  ierr = nf90_inq_dimid(ncid,'Nnode',dimid)
  call check2(ierr, 'inq_dimid Nnode')

  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%ncoord)
  call check2(ierr, 'inquire_dimension Node')

  if (nproc > 1) then
  ierr = nf90_inq_dimid(ncid,'mpi_nn',dimid)
  call check2(ierr, 'inq_dimid mpi_nn')

  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%mpi_nn)
  call check2(ierr, 'inquire_dimension mpi_nn')

  ierr = nf90_inq_dimid(ncid,'mpi_ne',dimid)
  call check2(ierr, 'inq_dimid mpi_ne')

  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%mpi_ne)
  call check2(ierr, 'inquire_dimension mpi_ne')

  ierr = nf90_inq_dimid(ncid,'mpi_nemax',dimid)
  call check2(ierr, 'inq_dimid mpi_nemax')

  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%mpi_nemax)
  call check2(ierr, 'inquire_dimension mpi_nemax')

  ierr = nf90_inq_dimid(ncid,'pinterfaces',dimid)
  call check2(ierr, 'inq_dimid pinterfaces')

  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%pinterfaces)
  call check2(ierr, 'inquire_dimension pinterfaces')
  end if

  allocate(this%coord    (     3,this%ncoord))
  allocate(this%elem     (     4,this%nelem))
  allocate(this%neigh    (Nfaces,this%nelem))
  allocate(this%face     (Nfaces,this%nelem))
  allocate(this%direction(Nfaces,this%nelem))
  allocate(this%bctype   (Nfaces,this%nelem))
  allocate(this%fluxtype   (Nfaces,this%nelem))
  allocate(this%elemtype (       this%nelem))
  ! media
  allocate(this%rho(this%nelem))
  allocate(this%vp (this%nelem))
  allocate(this%vs (this%nelem))
  allocate(this%zp (this%nelem))
  allocate(this%zs (this%nelem))
  allocate(this%mu (this%nelem))
  allocate(this%lam(this%nelem))
  allocate(this%vol(this%nelem))

  if (nproc > 1) then
  allocate(this%mpi_neighbor(this%mpi_nn))
  allocate(this%mpi_connection(this%mpi_nn,this%mpi_ne,2))
  allocate(this%mpi_ibool(Nfaces,this%mpi_nemax))
  allocate(this%mpi_interface(4,Nfaces,this%nelem))
  allocate(this%mpi_vp(this%pinterfaces,3))
  allocate(this%mpi_vs(this%pinterfaces,3))
  allocate(this%mpi_rho(this%pinterfaces,3))
  end if

  ierr = nf90_inq_varid(ncid,'node',varid)
  call check2(ierr, 'inq_varid node')

  ierr = nf90_get_var(ncid,varid,this%coord)
  call check2(ierr, 'get_var node')

  ierr = nf90_inq_varid(ncid,'elem',varid)
  call check2(ierr, 'inq_varid elem')

  ierr = nf90_get_var(ncid,varid,this%elem)
  call check2(ierr, 'get_var elem')

  ierr = nf90_inq_varid(ncid,'neighbor',varid)
  call check2(ierr, 'inq_varid neigh')

  ierr = nf90_get_var(ncid,varid,this%neigh)
  call check2(ierr, 'get_var neigh')

  ierr = nf90_inq_varid(ncid,'face',varid)
  call check2(ierr, 'inq_varid face')

  ierr = nf90_get_var(ncid,varid,this%face)
  call check2(ierr, 'get_var face')

  ierr = nf90_inq_varid(ncid,'direction',varid)
  call check2(ierr, 'inq_varid direction')

  ierr = nf90_get_var(ncid,varid,this%direction)
  call check2(ierr, 'get_var direction')

  ierr = nf90_inq_varid(ncid,'bctype',varid)
  call check2(ierr, 'inq_varid bctype')

  ierr = nf90_get_var(ncid,varid,this%bctype)
  call check2(ierr, 'get_var bctype')

  ierr = nf90_inq_varid(ncid,'fluxtype',varid)
  call check2(ierr, 'inq_varid fluxtype')

  ierr = nf90_get_var(ncid,varid,this%fluxtype)
  call check2(ierr, 'get_var fluxtype')

  ierr = nf90_inq_varid(ncid,'elemtype',varid)
  call check2(ierr, 'inq_varid elemtype')

  ierr = nf90_get_var(ncid,varid,this%elemtype)
  call check2(ierr, 'get_var elemtype')

  ierr = nf90_inq_varid(ncid,'rho',varid)
  call check2(ierr, 'inq_varid rho')

  ierr = nf90_get_var(ncid,varid,this%rho)
  call check2(ierr, 'get_var rho')

  ierr = nf90_inq_varid(ncid,'vp',varid)
  call check2(ierr, 'inq_varid vp')

  ierr = nf90_get_var(ncid,varid,this%vp)
  call check2(ierr, 'get_var vp')

  ierr = nf90_inq_varid(ncid,'vs',varid)
  call check2(ierr, 'inq_varid vs')

  ierr = nf90_get_var(ncid,varid,this%vs)
  call check2(ierr, 'get_var vs')

  if (nproc > 1) then
  ierr = nf90_inq_varid(ncid,'mpi_neighbor',varid)
  call check2(ierr, 'inq_varid mpi_neigh')

  ierr = nf90_get_var(ncid,varid,this%mpi_neighbor)
  call check2(ierr, 'get_var mpi_neigh')

  ierr = nf90_inq_varid(ncid,'mpi_connection',varid)
  call check2(ierr, 'inq_varid mpi_connection')

  ierr = nf90_get_var(ncid,varid,this%mpi_connection)
  call check2(ierr, 'get_var mpi_connection')

  ierr = nf90_inq_varid(ncid,'mpi_ibool',varid)
  call check2(ierr, 'inq_varid mpi_ibool')

  ierr = nf90_get_var(ncid,varid,this%mpi_ibool)
  call check2(ierr, 'get_var mpi_ibool')

  ierr = nf90_inq_varid(ncid,'mpi_interface',varid)
  call check2(ierr, 'inq_varid mpi_interface')

  ierr = nf90_get_var(ncid,varid,this%mpi_interface)
  call check2(ierr, 'get_var mpi_interface')

  ierr = nf90_inq_varid(ncid,'mpi_vp',varid)
  call check2(ierr, 'inq_varid mpi_vp')

  ierr = nf90_get_var(ncid,varid,this%mpi_vp)
  call check2(ierr, 'get_var mpi_vp')

  ierr = nf90_inq_varid(ncid,'mpi_vs',varid)
  call check2(ierr, 'inq_varid mpi_vs')

  ierr = nf90_get_var(ncid,varid,this%mpi_vs)
  call check2(ierr, 'get_var mpi_vs')

  ierr = nf90_inq_varid(ncid,'mpi_rho',varid)
  call check2(ierr, 'inq_varid mpi_rho')

  ierr = nf90_get_var(ncid,varid,this%mpi_rho)
  call check2(ierr, 'get_var mpi_rho')
  end if

  ierr = nf90_close(ncid)
  call check2(ierr, 'nf90_close '//trim(filename))

  !print*,'rank=',myrank,'ncoord=',this%ncoord,'nelem=',this%nelem!,'fault face=',this%nfault_face
  !print*,'rank=',myrank,'vp=',minval(this%vp),'-',maxval(this%vp)
  !print*,'rank=',myrank,'vs=',minval(this%vs),'-',maxval(this%vs)
  !print*,'rank=',myrank,'rho=',minval(this%rho),'-',maxval(this%rho)
  !stop 2

  ! write vtk data that can be read by paraview
  if (export_media .eq. 1) then
    write(filename,'(a,a,i6.6,a)') trim(mesh_dir),'/rho',myrank,'.vtk'
    write(*,*) "write ", filename
    call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%rho)

    write(filename,'(a,a,i6.6,a)') trim(mesh_dir),'/vp',myrank,'.vtk'
    write(*,*) "write ", filename
    call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%vp)

    write(filename,'(a,a,i6.6,a)') trim(mesh_dir),'/vs',myrank,'.vtk'
    write(*,*) "write ", filename
    call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%vs)
  end if

end subroutine

!@subroutine read_mesh(mesh)
!@  implicit none
!@  type(meshvar),intent(inout) :: mesh
!@  integer :: i,j!,k
!@  character(len=128) :: basename = 'mesh1'
!@
!@  ! read in nodes
!@  open(100,file=trim(basename)//'.node')
!@  read(100,*) mesh%ncoord,i
!@  allocate(mesh%coord(3,mesh%ncoord))
!@  do j = 1,mesh%ncoord
!@    read(100,*) (mesh%coord(i,j),i=1,3)
!@  end do
!@  close(100)
!@
!@  ! read in elements
!@  open(100,file=trim(basename)//'.elem')
!@  read(100,*) mesh%nelem,i
!@  allocate(mesh%elem(4,mesh%nelem))
!@  do j = 1,mesh%nelem
!@    read(100,*) (mesh%elem(i,j),i=1,4)
!@  enddo
!@  close(100)
!@
!@  ! read in neighbours
!@  open(100,file=trim(basename)//'.neigh')
!@  read(100,*) i,j
!@  allocate(mesh%neigh(Nfaces,mesh%nelem))
!@  do j = 1,mesh%nelem
!@    read(100,*) (mesh%neigh(i,j),i=1,4)
!@  enddo
!@  close(100)
!@
!@  ! read in face number of neighbour elements
!@  open(100,file=trim(basename)//'.face')
!@  read(100,*) i,j
!@  allocate(mesh%face(Nfaces,mesh%nelem))
!@  do j = 1,mesh%nelem
!@    read(100,*) (mesh%face(i,j),i=1,4)
!@  enddo
!@  close(100)
!@
!@  ! read in face reordering
!@  open(100,file=trim(basename)//'.direction')
!@  read(100,*) i,j
!@  allocate(mesh%direction(Nfaces,mesh%nelem))
!@  do j = 1,mesh%nelem
!@    read(100,*) (mesh%direction(i,j),i=1,4)
!@  enddo
!@  close(100)
!@
!@  allocate(mesh%bctype(Nfaces,mesh%nelem))
!@
!@  if(.false.) then
!@    print*,'ncoord=',mesh%ncoord
!@    do j = 1,mesh%ncoord
!@      print*,(mesh%coord(i,j),i=1,3)
!@    end do
!@
!@    print*,'nelem=',mesh%nelem
!@    do j = 1,mesh%nelem
!@      print*,(mesh%elem(i,j),i=1,4)
!@    enddo
!@
!@    print*,'neigh=',mesh%nelem
!@    do j = 1,mesh%nelem
!@      print*,(mesh%neigh(i,j),i=1,4)
!@    enddo
!@
!@    print*,'face=',mesh%nelem
!@    do j = 1,mesh%nelem
!@      print*,(mesh%face(i,j),i=1,4)
!@    enddo
!@
!@    print*,'direction=',mesh%nelem
!@    do j = 1,mesh%nelem
!@      print*,(mesh%direction(i,j),i=1,4)
!@    enddo
!@  end if
!@
!@
!@end subroutine

end module
