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

program main

use netcdf

implicit none

integer,parameter :: RKIND = 8
integer,parameter :: Nfaces = 4
integer,parameter :: nsize = 100 ! elements share one node
integer,parameter :: max_neighbor = 12 ! max neighbors per element
integer,parameter :: BC_FAULT = 100

type :: meshvar
  ! basic parameters
  integer :: nelem,ncoord,nfault_elem
  integer :: nglob

  ! mpi parameters
  integer :: mpi_nn
  integer :: mpi_ne
  integer :: mpi_nnmax
  integer :: mpi_nemax
  integer :: pinterfaces

  ! basic connections
  real(kind=RKIND),dimension(:,:),pointer :: coord=>null()
  integer,dimension(:,:),pointer :: elem => null()
  integer,dimension(:,:),pointer :: neighbor => null()
  integer,dimension(:,:),pointer :: face => null()
  integer,dimension(:,:),pointer :: direction => null()
  integer,dimension(:,:),pointer :: bctype => null()
  integer,dimension(:,:),pointer :: fluxtype => null()
  integer,dimension(:),pointer :: elemtype => null()

  ! wave to fault elem index
  integer,dimension(:),pointer :: wave2fault => null()
  integer,dimension(:),pointer :: fault2wave => null()

  ! metis partition
  integer,dimension(:),pointer :: part => null()

  ! mpi
  integer,dimension(:),pointer :: mpi_neighbor => null()
  integer,dimension(:,:,:),pointer :: mpi_connection => null()
  integer,dimension(:,:),pointer :: mpi_ibool => null()
  integer,dimension(:,:,:),pointer :: mpi_interface => null()
  integer,dimension(:),pointer :: mpi_ninterface => null()
  real(kind=rkind),dimension(:,:),pointer :: mpi_vp => null()
  real(kind=rkind),dimension(:,:),pointer :: mpi_vs => null()
  real(kind=rkind),dimension(:,:),pointer :: mpi_rho => null()

  ! media
  real(kind=rkind),dimension(:),pointer :: vp => null()
  real(kind=rkind),dimension(:),pointer :: vs => null()
  real(kind=rkind),dimension(:),pointer :: rho => null()

  ! initial stress
  real(kind=rkind),dimension(:,:,:),pointer :: Tx0 => null()
  real(kind=rkind),dimension(:,:,:),pointer :: Ty0 => null()
  real(kind=rkind),dimension(:,:,:),pointer :: Tz0 => null()
  real(kind=rkind),dimension(:,:,:),pointer :: dTx0 => null()
  real(kind=rkind),dimension(:,:,:),pointer :: dTy0 => null()
  real(kind=rkind),dimension(:,:,:),pointer :: dTz0 => null()

  ! friction database
  real(kind=rkind),dimension(:,:,:),pointer :: mu_s => null()
  real(kind=rkind),dimension(:,:,:),pointer :: mu_d => null()
  real(kind=rkind),dimension(:,:,:),pointer :: Dc => null()
  real(kind=rkind),dimension(:,:,:),pointer :: C0 => null()
  ! rate state
  real(kind=rkind),dimension(:,:,:),pointer :: a => null()
  real(kind=rkind),dimension(:,:,:),pointer :: b => null()
  real(kind=rkind),dimension(:,:,:),pointer :: Vw => null()
  real(kind=rkind),dimension(:,:,:),pointer :: state => null()
  ! thermal pressurization
  real(kind=rkind),dimension(:,:,:),pointer :: TP_hy => null()


end type

type paras
  integer :: nproc
end type


integer, dimension(:),allocatable :: part ! contains to which rank an element belongs
integer, dimension(:), allocatable :: glob2loc_elmnts ! translates global to local rank numbering of elements
integer, dimension(:,:), allocatable :: loc2glob_elemnts ! translates local numbering back to global numbering
integer, dimension(:), allocatable :: num_loc,num_loc_fault
integer, dimension(:), allocatable :: nodes_elmnts, nnodes_elmnts ! global numbering of nodes, number of elements per node
integer, dimension(:), allocatable :: glob2loc_nodes_nparts ! number of a node in overall sum _size_glob2loc_nodes_
integer, dimension(:), allocatable :: glob2loc_nodes_parts ! rank where node of overall numbering is located
integer, dimension(:), allocatable :: glob2loc_nodes,glob2loc_nodes2  ! local numbering of nodes in each rank depending on global overall numbering
integer :: size_glob2loc_nodes ! counter for overall sum of nodes if nodes belonging to multiple ranks counted multiple
integer, dimension(:), allocatable :: part_nodes, num_parts ! shows if a node belongs to a rank, number of nodes in a rank
integer, dimension(4) :: loc_nodes ! local number of a node in its rank
integer, dimension(:,:), allocatable :: icom ! matrix with entry 1 if two ranks connected, 0 otherwise
integer, dimension(:), allocatable :: elmnts ! vector containing nodes of all elements
integer, dimension(:), allocatable :: tempv

type(paras) :: par
type(meshvar) :: this
type(meshVar), dimension(:), allocatable :: db ! container for mesh variables of every rank for MPI run

integer :: i,j,k,l,m,c,is,ie,ief,lf,in,idb,iz,ee,iproc,temp1,temp2
integer :: num_glob,num_part
logical :: isfault
character(len=128) :: filename
character(len=128) :: arg

call get_command_argument(1, arg)
if(LEN_TRIM(arg) == 0) then
  print*,'mesh filename required!'
  print*,'Usage: ./part_mesh mesh.nc'
  STOP 2
end if

print*,trim(arg)

!call read_mesh(this,'mesh.nc')
call read_mesh(this,trim(arg))
allocate(part(this%nelem))
do i = 1,this%nelem
part(i) = this%part(i)
end do
par%nproc = maxval(part)
print*,'nproc=',par%nproc
!print*,this%part

allocate(db(par%nproc))

allocate(elmnts(this%nelem*4))
k=1
do i=1,this%nelem
  do j=1,4
    elmnts(k)=this%elem(j,i)
    k=k+1
  end do
end do
! create glob2loc_elem
! inspired by the specfem partition (SPECFEM3D 2.0 in subroutine part_decompose_mesh_SCOTCH.f90
! written by Komatitsch at al. published under GPL taken from www.geodynamics.org)
! be careful, metis gives parts begining from 0 or from 1 depending on the version,
!   so I compile a local version, here 4.0.3
allocate(glob2loc_elmnts(this%nelem))
glob2loc_elmnts(:) = 0

allocate(num_loc(par%nproc))
allocate(num_loc_fault(par%nproc))
num_loc(:) = 1
num_loc_fault(:) = 1

! build a vector with the local numbering of the elements
do num_glob=1,this%nelem
  num_part=part(num_glob)
  glob2loc_elmnts(num_glob) = num_loc(num_part)
  num_loc(num_part) = num_loc(num_part) + 1
  if (this%wave2fault(num_glob)>0) then
    num_loc_fault(num_part)=num_loc_fault(num_part)+1
  end if
end do
do iproc = 1,par%nproc
  print*,'iproc=',iproc,'num_loc=',num_loc(iproc),'num_loc_fault=',num_loc_fault(iproc)-1
end do

! create local node numbering
! 2 vectors, nnodes with the positions of the elements in the vector nodes_elements, similar to the metis numbering
allocate(nnodes_elmnts(this%ncoord))
allocate(nodes_elmnts(this%ncoord*nsize))
nnodes_elmnts(:)=0
nodes_elmnts(:)=0
!print*,'elmnts=',elmnts
do i=1, 4*this%nelem
  ! nodes is like a matrix with nodes as rows and nsize elements as colums
  nodes_elmnts((elmnts(i)-1)*nsize + nnodes_elmnts(elmnts(i))+1) = 1+(i-1)/4
  nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) +1
end do

print*,'create local node numbering ...'
! create the local node numbering
allocate(glob2loc_nodes_nparts(this%ncoord+1))
allocate(part_nodes(par%nproc), num_parts(par%nproc))

size_glob2loc_nodes = 1
part_nodes(:) = 0

do in=1,this%ncoord
  glob2loc_nodes_nparts(in) = size_glob2loc_nodes
  do ie = 1, nnodes_elmnts(in)
    ! shows to which partitions a node belongs
    part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1
  end do

  do num_part=1,par%nproc
    if (part_nodes(num_part) == 1) then
      ! get number of nodes in the global array, there might be some nodes which count double at the interfaces
      size_glob2loc_nodes = size_glob2loc_nodes +1
      part_nodes(num_part) = 0
    end if
  end do
  glob2loc_nodes_nparts(this%ncoord+1) = size_glob2loc_nodes
end do

allocate(glob2loc_nodes_parts(glob2loc_nodes_nparts(this%ncoord+1)-1))
allocate(glob2loc_nodes(glob2loc_nodes_nparts(this%ncoord+1)-1))
allocate(glob2loc_nodes2(this%ncoord))

!glob2loc_nodes_parts(:) = -999

glob2loc_nodes(:) = 1
part_nodes(:) = 0
num_parts(:)=1
size_glob2loc_nodes = 1

do in=1,this%ncoord
  do ie = 1, nnodes_elmnts(in)
    part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1
  end do

  do num_part=1,par%nproc
    if (part_nodes(num_part) == 1) then
      ! build arrays with local nodes ordering
      glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
      glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
      size_glob2loc_nodes = size_glob2loc_nodes +1
      num_parts(num_part) = num_parts(num_part) +1
      part_nodes(num_part) = 0
    end if
  end do
end do

print*,'build databases ...'
!"--------------------------------------------------------------------------------------"
! build Databases
!"--------------------------------------------------------------------------------------"
allocate(loc2glob_elemnts(par%nproc,maxval(num_loc(:)-1)))
do iproc=1,par%nproc
  if(mod(iproc,max(par%nproc/20,1))==0) print*,iproc,'/',par%nproc
  ! build local to global numbering
  do ie=1,this%nelem
    if (part(ie)== iproc) then
      loc2glob_elemnts(iproc,glob2loc_elmnts(ie))=ie
    endif
  enddo

  !prepare database for every partition
  !db(iproc)%nglob = (num_loc(iproc)-1)*Np
  db(iproc)%nelem = num_loc(iproc)-1
  db(iproc)%ncoord = num_parts(iproc)-1
  db(iproc)%nfault_elem = num_loc_fault(iproc)-1
  allocate(db(iproc)%coord(3,db(iproc)%ncoord))
  allocate(db(iproc)%elem(4,db(iproc)%nelem))
  allocate(db(iproc)%neighbor(4,db(iproc)%nelem))
  allocate(db(iproc)%face(4,db(iproc)%nelem))
  allocate(db(iproc)%direction(4,db(iproc)%nelem))
  allocate(db(iproc)%bctype(4,db(iproc)%nelem))
  allocate(db(iproc)%fluxtype(4,db(iproc)%nelem))
  allocate(db(iproc)%elemtype(db(iproc)%nelem))
  allocate(db(iproc)%vp(db(iproc)%nelem))
  allocate(db(iproc)%vs(db(iproc)%nelem))
  allocate(db(iproc)%rho(db(iproc)%nelem))
  allocate(db(iproc)%mpi_interface(4,4,db(iproc)%nelem))

  !db(1)%nfault_elem=0
  allocate(db(iproc)%Tx0 (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%Ty0 (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%Tz0 (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%dTx0 (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%dTy0 (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%dTz0 (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%mu_s(3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%mu_d(3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%Dc  (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%C0  (3,Nfaces,db(iproc)%nfault_elem))
  ! rate state
  allocate(db(iproc)%a   (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%b   (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%Vw  (3,Nfaces,db(iproc)%nfault_elem))
  allocate(db(iproc)%state(3,Nfaces,db(iproc)%nfault_elem))
  ! thermal pressurization
  allocate(db(iproc)%TP_hy(3,Nfaces,db(iproc)%nfault_elem))

  ! build local arrays and coordinates
  do i=1, this%ncoord
    do j= glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
      if (glob2loc_nodes_parts(j) == iproc) then
        k=glob2loc_nodes(j)
        db(iproc)%coord(1,k)=this%coord(1,i)
        db(iproc)%coord(2,k)=this%coord(2,i)
        db(iproc)%coord(3,k)=this%coord(3,i)
        !db(iproc)%loc2glob_nodes(k)=i
      end if
    end do
  end do

  ! build local elements
  do i=1, this%nelem
    if (part(i) == iproc) then
      do j=1,4
        l=elmnts((i-1)*4+j)
        do k=glob2loc_nodes_nparts(l), glob2loc_nodes_nparts(l+1)-1
          if (glob2loc_nodes_parts(k) == iproc) then
            loc_nodes(j) = glob2loc_nodes(k)
          end if
        end do
      end do
      k = glob2loc_elmnts(i)
      db(iproc)%elem(1,k) = loc_nodes(1)
      db(iproc)%elem(2,k) = loc_nodes(2)
      db(iproc)%elem(3,k) = loc_nodes(3)
      db(iproc)%elem(4,k) = loc_nodes(4)
    end if
  end do

  !"---------------------------------------------------"
  ! find mpi neighbors and the neighbor element
  !"---------------------------------------------------"

  db(iproc)%neighbor(:,:) = 0
  db(iproc)%mpi_interface(:,:,:) = 0
  db(iproc)%pinterfaces=0
  do ie=1,this%nelem
    do j=1,4
      k=this%neighbor(j,ie)
      if (k>0) then
        if (part(ie) == iproc) then
          if (part(k) == iproc) then
            l = glob2loc_elmnts(ie)
            m = glob2loc_elmnts(k)
            db(iproc)%neighbor(j,l) = m
          else
            l = glob2loc_elmnts(ie)
            m = glob2loc_elmnts(k)
            db(iproc)%neighbor(j,l) = -1 ! means mpi neighbor
            db(iproc)%mpi_interface(1,j,l) = part(k) ! neighbor partition
            db(iproc)%mpi_interface(2,j,l) = m ! element in partition
            db(iproc)%pinterfaces=db(iproc)%pinterfaces+1
          end if
        end if
      end if
    end do
  end do

  ! add physical parameters to the databases
  do ie=1, this%nelem
    if (part(ie)== iproc) then
      l = glob2loc_elmnts(ie)
      db(iproc)%face(:,l) = this%face(:,ie)
      db(iproc)%direction(:,l) = this%direction(:,ie)
      db(iproc)%bctype(:,l) = this%bctype(:,ie)
      db(iproc)%fluxtype(:,l) = this%fluxtype(:,ie)
      db(iproc)%elemtype(l) = this%elemtype(ie)
      db(iproc)%rho(l)=this%rho(ie)
      db(iproc)%vp(l)=this%vp(ie)
      db(iproc)%vs(l)=this%vs(ie)
    end if
  end do

  allocate(db(iproc)%fault2wave(db(iproc)%nfault_elem))
  allocate(db(iproc)%wave2fault(db(iproc)%nelem))

  db(iproc)%wave2fault(:) = 0

  k = 0
  do ie = 1,db(iproc)%nelem
    isfault=.false.
    do is = 1,Nfaces
      if (db(iproc)%bctype(is,ie) >= BC_FAULT) then
        isfault=.true.
      end if
    end do
    if (isfault) then
      k = k + 1
      db(iproc)%fault2wave(k) = ie
      db(iproc)%wave2fault(ie) = k
    end if
  end do

  ! add physical parameters to the databases
  do ief=1, this%nfault_elem
    ie = this%fault2wave(ief)
    if (part(ie)== iproc) then
      l = glob2loc_elmnts(ie)
      lf = db(iproc)%wave2fault(l)
      db(iproc)%Tx0(:,:,lf)=this%Tx0(:,:,ief)
      db(iproc)%Ty0(:,:,lf)=this%Ty0(:,:,ief)
      db(iproc)%Tz0(:,:,lf)=this%Tz0(:,:,ief)
      db(iproc)%dTx0(:,:,lf)=this%dTx0(:,:,ief)
      db(iproc)%dTy0(:,:,lf)=this%dTy0(:,:,ief)
      db(iproc)%dTz0(:,:,lf)=this%dTz0(:,:,ief)
      db(iproc)%mu_s(:,:,lf)=this%mu_s(:,:,ief)
      db(iproc)%mu_d(:,:,lf)=this%mu_d(:,:,ief)
      db(iproc)%Dc(:,:,lf)=this%Dc(:,:,ief)
      db(iproc)%C0(:,:,lf)=this%C0(:,:,ief)
      ! rate state
      db(iproc)%a(:,:,lf)=this%a(:,:,ief)
      db(iproc)%b(:,:,lf)=this%b(:,:,ief)
      db(iproc)%Vw(:,:,lf)=this%Vw(:,:,ief)
      db(iproc)%state(:,:,lf)=this%state(:,:,ief)
      ! thermal pressurization
      db(iproc)%TP_hy(:,:,lf)=this%TP_hy(:,:,ief)
    end if
  end do


end do ! nproc

print*,'create mpi interfaces ...'
!"--------------------------------------------------------------------------------------"
! create mpi interfaces
!"--------------------------------------------------------------------------------------"
allocate(icom(par%nproc,par%nproc))
icom(:,:) = 0
do iproc=1,par%nproc
  do ie=1,db(iproc)%nelem
    do i=1,4
      if (db(iproc)%mpi_interface(1,i,ie) > 0) then
        icom(iproc,db(iproc)%mpi_interface(1,i,ie))=1
      end if
    end do
  end do
  db(iproc)%mpi_nn=sum(icom(iproc,:))
end do

print*, 'set up mpi_neighbor array'
do iproc=1,par%nproc
  allocate(db(iproc)%mpi_neighbor(db(iproc)%mpi_nn))
  c=1
  do i=1,par%nproc
    if (icom(iproc,i) == 1) then
      db(iproc)%mpi_neighbor(c)=i
      c=c+1
    end if
  end do
end do

! how many mpi interfaces do we have?
do iproc=1,par%nproc
  allocate(tempv(db(iproc)%mpi_nn))
  allocate(db(iproc)%mpi_ninterface(db(iproc)%mpi_nn))
  tempv(:)=0
  do i=1,db(iproc)%mpi_nn
    in=db(iproc)%mpi_neighbor(i)
    do ie=1,db(iproc)%nelem
      do is=1,4
        if ( (db(iproc)%neighbor(is,ie) == -1) .and. (db(iproc)%mpi_interface(1,is,ie) == in)) then
          tempv(i)=tempv(i)+1
        end if
      end do
    end do
!   write(*,*) "proc",iproc," found", tempv(i), "elements on interface between", iproc ,"and", in
    db(iproc)%mpi_ninterface(i)=tempv(i)
  end do
  deallocate(tempv)
end do ! nproc

temp2=0
do iproc=1,par%nproc
  temp1=db(iproc)%nelem
  if (temp2<temp1) temp2=temp1
end do

! set up mpi_ibool with max mpi neighbors
do iproc=1,par%nproc
  db(iproc)%mpi_nemax=temp2
  allocate(db(iproc)%mpi_ibool(4,db(iproc)%mpi_nemax))
  db(iproc)%mpi_ibool(:,:)=0
end do

temp2=0
do iproc=1,par%nproc
  temp1=maxval(db(iproc)%mpi_ninterface(:))
  if (temp2<temp1) temp2=temp1
end do
!"--------------------------------------------------------------------------------------"
! create main mpi reconnection arrays.
!"--------------------------------------------------------------------------------------"
do iproc=1,par%nproc
  db(iproc)%mpi_nnmax=maxval(db(:)%mpi_nn)
  !allocate(db(iproc)%mpi_icon(db(iproc)%mpi_nnmax))
  !db(iproc)%mpi_icon(:) = 0
  db(iproc)%mpi_ne=temp2
  allocate(db(iproc)%mpi_connection(db(iproc)%mpi_nn,db(iproc)%mpi_ne,2))
  db(iproc)%mpi_connection(:,:,:)=0
  do i=1,db(iproc)%mpi_nn
    in=db(iproc)%mpi_neighbor(i)
    c=1
    do ie=1,db(iproc)%nelem
      do is=1,4
        if ( (db(iproc)%neighbor(is,ie) == -1) .and. (db(iproc)%mpi_interface(1,is,ie) == in)) then
          db(iproc)%mpi_connection(i,c,1)=ie ! mpi interface to global element
          db(iproc)%mpi_connection(i,c,2)=is ! which side of the element is the interface
          db(iproc)%mpi_interface(3,is,ie)=c ! global element to local mpi interface
          db(iproc)%mpi_interface(4,is,ie)=i ! which local interface?
          c=c+1
        end if
      end do
    end do
  end do
end do ! nproc

!create Impedance mpi connection array
do iproc = 1, par%nproc
  allocate(db(iproc)%mpi_vp(db(iproc)%pinterfaces, 3))
  allocate(db(iproc)%mpi_vs(db(iproc)%pinterfaces, 3))
  allocate(db(iproc)%mpi_rho(db(iproc)%pinterfaces, 3))
  i = 1
  do ie = 1, db(iproc)%nelem
    do is = 1, Nfaces
      if (db(iproc)%neighbor(is,ie) == -1) then
        idb = db(iproc)%mpi_interface(1,is,ie)
        iz  = db(iproc)%mpi_interface(2,is,ie)
        db(iproc)%mpi_vp(i,1) = db(idb)%vp(iz)
        db(iproc)%mpi_vp(i,2) = idb
        db(iproc)%mpi_vp(i,3) = iz
        db(iproc)%mpi_vs(i,1) = db(idb)%vs(iz)
        db(iproc)%mpi_vs(i,2) = idb
        db(iproc)%mpi_vs(i,3) = iz
        db(iproc)%mpi_rho(i,1) = db(idb)%rho(iz)
        db(iproc)%mpi_rho(i,2) = idb
        db(iproc)%mpi_rho(i,3) = iz
        i = i+1
      end if
    end do
  end do
end do

do iproc=1,par%nproc
  do i=1,db(iproc)%mpi_nn
    l=db(iproc)%mpi_neighbor(i)
    do ie=1,db(iproc)%mpi_ne
      if ( db(iproc)%mpi_connection(i,ie,1) >0) then
        is=db(iproc)%mpi_connection(i,ie,2)
        ee=db(iproc)%mpi_connection(i,ie,1)
        k=db(iproc)%mpi_interface(2,is,ee)
        c = db(l)%mpi_interface(3,db(iproc)%face(is,ee),k)
        !db(iproc)%mpi_icon(i)=db(l)%mpi_interface(4,db(iproc)%face(is,ee),k)
        ! correct mpi boundarys
        db(iproc)%mpi_ibool(is,ee)=c
      end if
    end do
  end do
end do

! write databases
do iproc=1,par%nproc
  write(filename,"('data/meshVar',i6.6,'.nc')") iproc-1
  print*,trim(filename)
  call write_mesh(db(iproc),trim(filename),par%nproc)
end do

contains

!subroutine check(status)
!  implicit none
!  integer, intent ( in) :: status
!
!  if(status /= nf90_noerr) then
!    print *, trim(nf90_strerror(status))!,' in file ',__FILE__,' line ',__LINE__
!    stop "Stopped"
!  end if
!end subroutine

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

subroutine read_mesh(this,filename)
  implicit none
  type(meshvar) :: this
  character(len=*) :: filename
  character(len=128) :: namestr
  !integer :: ios
  integer :: ierr,ncid,dimid,varid
  !integer :: Nelem,Nnode

  !write(filename,'(a,i6.6,a)') 'data/meshVar',myrank,'.nc'
  ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid)
  call check2(ierr,'nf90_open read_mesh')

  ierr = nf90_inq_dimid(ncid,'Nelem',dimid)
  call check2(ierr,'inq_dimid Nelem')
  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%nelem)
  call check2(ierr,'inq_dimension Nelem')

  ierr = nf90_inq_dimid(ncid,'Nnode',dimid)
  call check2(ierr,'inq_dimid Nnode')
  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%ncoord)
  call check2(ierr,'inq_dimension Nnode')

  ierr = nf90_inq_dimid(ncid,'nfault_elem',dimid)
  call check2(ierr,'inq_dimid nfault_elem')
  ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%nfault_elem)
  call check2(ierr,'inq_dimension nfault_elem')

  !ierr = nf90_inq_dimid(ncid,'mpi_nn',dimid)
  !if (ierr /= nf90_noerr) stop 16
  !ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%mpi_nn)
  !if (ierr /= nf90_noerr) stop 17

  !ierr = nf90_inq_dimid(ncid,'mpi_ne',dimid)
  !if (ierr /= nf90_noerr) stop 18
  !ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%mpi_ne)
  !if (ierr /= nf90_noerr) stop 19

  !ierr = nf90_inq_dimid(ncid,'mpi_nemax',dimid)
  !if (ierr /= nf90_noerr) stop 20
  !ierr = nf90_inquire_dimension(ncid,dimid,namestr,this%mpi_nemax)
  !if (ierr /= nf90_noerr) stop 21

  allocate(this%coord    (     3,this%ncoord))
  allocate(this%elem     (     4,this%nelem))
  allocate(this%neighbor (Nfaces,this%nelem))
  allocate(this%face     (Nfaces,this%nelem))
  allocate(this%direction(Nfaces,this%nelem))
  allocate(this%bctype   (Nfaces,this%nelem))
  allocate(this%fluxtype   (Nfaces,this%nelem))
  allocate(this%elemtype (       this%nelem))
  allocate(this%vp       (       this%nelem))
  allocate(this%vs       (       this%nelem))
  allocate(this%rho      (       this%nelem))

  allocate(this%part     (       this%nelem))

  allocate(this%Tx0 (3,4,this%nfault_elem))
  allocate(this%Ty0 (3,4,this%nfault_elem))
  allocate(this%Tz0 (3,4,this%nfault_elem))
  allocate(this%dTx0 (3,4,this%nfault_elem))
  allocate(this%dTy0 (3,4,this%nfault_elem))
  allocate(this%dTz0 (3,4,this%nfault_elem))
  allocate(this%mu_s(3,4,this%nfault_elem))
  allocate(this%mu_d(3,4,this%nfault_elem))
  allocate(this%Dc  (3,4,this%nfault_elem))
  allocate(this%C0  (3,4,this%nfault_elem))
  ! rate state
  allocate(this%a   (3,4,this%nfault_elem))
  allocate(this%b   (3,4,this%nfault_elem))
  allocate(this%Vw  (3,4,this%nfault_elem))
  allocate(this%state(3,4,this%nfault_elem))
  ! thermal pressurization
  allocate(this%TP_hy(3,4,this%nfault_elem))

  allocate(this%wave2fault(this%nelem))
  allocate(this%fault2wave(this%nfault_elem))


  !allocate(this%mpi_neighbor(this%mpi_nn))
  !allocate(this%mpi_connection(this%mpi_nn,this%mpi_ne,2))
  !allocate(this%mpi_ibool(Nfaces,this%mpi_nemax))
  !allocate(this%mpi_interface(4,Nfaces,this%nelem))

  ierr = nf90_inq_varid(ncid,'node',varid)
  call check2(ierr,'inq_varid node')
  ierr = nf90_get_var(ncid,varid,this%coord)
  call check2(ierr,'get_var node')

  ierr = nf90_inq_varid(ncid,'elem',varid)
  call check2(ierr,'inq_varid elem')
  ierr = nf90_get_var(ncid,varid,this%elem)
  call check2(ierr,'get_var elem')

  ierr = nf90_inq_varid(ncid,'neighbor',varid)
  call check2(ierr,'inq_varid neighbor')
  ierr = nf90_get_var(ncid,varid,this%neighbor)
  call check2(ierr,'get_var neighbor')

  ierr = nf90_inq_varid(ncid,'face',varid)
  call check2(ierr,'inq_varid face')
  ierr = nf90_get_var(ncid,varid,this%face)
  call check2(ierr,'get_var face')

  ierr = nf90_inq_varid(ncid,'direction',varid)
  call check2(ierr,'inq_varid direction')
  ierr = nf90_get_var(ncid,varid,this%direction)
  call check2(ierr,'get_var direction')

  ierr = nf90_inq_varid(ncid,'bctype',varid)
  call check2(ierr,'inq_varid bctype')
  ierr = nf90_get_var(ncid,varid,this%bctype)
  call check2(ierr,'get_var bctype')

  ierr = nf90_inq_varid(ncid,'fluxtype',varid)
  call check2(ierr,'inq_varid fluxtype')
  ierr = nf90_get_var(ncid,varid,this%fluxtype)
  call check2(ierr,'get_var fluxtype')

  ierr = nf90_inq_varid(ncid,'elemtype',varid)
  call check2(ierr,'inq_varid elemtype')
  ierr = nf90_get_var(ncid,varid,this%elemtype)
  call check2(ierr,'get_var elemtype')

  ierr = nf90_inq_varid(ncid,'vp',varid)
  call check2(ierr,'inq_varid vp')
  ierr = nf90_get_var(ncid,varid,this%vp)
  call check2(ierr,'get_var vp')

  ierr = nf90_inq_varid(ncid,'vs',varid)
  call check2(ierr,'inq_varid vs')
  ierr = nf90_get_var(ncid,varid,this%vs)
  call check2(ierr,'get_var vs')

  ierr = nf90_inq_varid(ncid,'rho',varid)
  call check2(ierr,'inq_varid rho')
  ierr = nf90_get_var(ncid,varid,this%rho)
  call check2(ierr,'get_var rho')

  ierr = nf90_inq_varid(ncid,'part',varid)
  call check2(ierr,'inq_varid part')
  ierr = nf90_get_var(ncid,varid,this%part)
  call check2(ierr,'get_var part')

  !ierr = nf90_inq_varid(ncid,'mpi_neighbor',varid)
  !ierr = nf90_get_var(ncid,varid,this%mpi_neighbor)

  !ierr = nf90_inq_varid(ncid,'mpi_connection',varid)
  !ierr = nf90_get_var(ncid,varid,this%mpi_connection)

  !ierr = nf90_inq_varid(ncid,'mpi_ibool',varid)
  !ierr = nf90_get_var(ncid,varid,this%mpi_ibool)

  !ierr = nf90_inq_varid(ncid,'mpi_interface',varid)
  !ierr = nf90_get_var(ncid,varid,this%mpi_interface)

  ierr = nf90_inq_varid(ncid,'Tx0',varid)
  call check2(ierr,'inq_varid Tx0')
  ierr = nf90_get_var(ncid,varid,this%Tx0)
  call check2(ierr,'get_var Tx0')

  ierr = nf90_inq_varid(ncid,'Ty0',varid)
  call check2(ierr,'inq_varid Ty0')
  ierr = nf90_get_var(ncid,varid,this%Ty0)
  call check2(ierr,'get_var Ty0')

  ierr = nf90_inq_varid(ncid,'Tz0',varid)
  call check2(ierr,'inq_varid Tz0')
  ierr = nf90_get_var(ncid,varid,this%Tz0)
  call check2(ierr,'get_var Tz0')

  ierr = nf90_inq_varid(ncid,'dTx0',varid)
  call check2(ierr,'inq_varid dTx0')
  ierr = nf90_get_var(ncid,varid,this%dTx0)
  call check2(ierr,'get_var dTx0')

  ierr = nf90_inq_varid(ncid,'dTy0',varid)
  call check2(ierr,'inq_varid dTy0')
  ierr = nf90_get_var(ncid,varid,this%dTy0)
  call check2(ierr,'get_var dTy0')

  ierr = nf90_inq_varid(ncid,'dTz0',varid)
  call check2(ierr,'inq_varid dTz0')
  ierr = nf90_get_var(ncid,varid,this%dTz0)
  call check2(ierr,'get_var dTz0')

  ierr = nf90_inq_varid(ncid,'mu_s',varid)
  call check2(ierr,'inq_varid mu_s')
  ierr = nf90_get_var(ncid,varid,this%mu_s)
  call check2(ierr,'get_var mu_s')

  ierr = nf90_inq_varid(ncid,'mu_d',varid)
  call check2(ierr,'inq_varid mu_d')
  ierr = nf90_get_var(ncid,varid,this%mu_d)
  call check2(ierr,'get_var mu_d')

  ierr = nf90_inq_varid(ncid,'Dc',varid)
  call check2(ierr,'inq_varid Dc')
  ierr = nf90_get_var(ncid,varid,this%Dc)
  call check2(ierr,'get_var Dc')

  ierr = nf90_inq_varid(ncid,'C0',varid)
  call check2(ierr,'inq_varid C0')
  ierr = nf90_get_var(ncid,varid,this%C0)
  call check2(ierr,'get_var C0')

  ! rate state
  ierr = nf90_inq_varid(ncid,'a',varid)
  call check2(ierr,'inq_varid a')
  ierr = nf90_get_var(ncid,varid,this%a)
  call check2(ierr,'get_var a')

  ierr = nf90_inq_varid(ncid,'b',varid)
  call check2(ierr,'inq_varid b')
  ierr = nf90_get_var(ncid,varid,this%b)
  call check2(ierr,'get_var b')

  ierr = nf90_inq_varid(ncid,'Vw',varid)
  call check2(ierr,'inq_varid Vw')
  ierr = nf90_get_var(ncid,varid,this%Vw)
  call check2(ierr,'get_var Vw')

  ierr = nf90_inq_varid(ncid,'state',varid)
  call check2(ierr,'inq_varid state')
  ierr = nf90_get_var(ncid,varid,this%state)
  call check2(ierr,'get_var state')

  ! thermal pressurization
  ierr = nf90_inq_varid(ncid,'TP_hy',varid)
  call check2(ierr,'inq_varid TP_hy')
  ierr = nf90_get_var(ncid,varid,this%TP_hy)
  call check2(ierr,'get_var TP_hy')

  ierr = nf90_inq_varid(ncid,'wave2fault',varid)
  call check2(ierr,'inq_varid wave2fault')
  ierr = nf90_get_var(ncid,varid,this%wave2fault)
  call check2(ierr,'get_var wave2fault')

  ierr = nf90_inq_varid(ncid,'fault2wave',varid)
  call check2(ierr,'inq_varid fault2wave')
  ierr = nf90_get_var(ncid,varid,this%fault2wave)
  call check2(ierr,'get_var fault2wave')

  !================================================
  ierr = nf90_close(ncid)
  !================================================
  print*,&
  'Nelem=',this%Nelem,&
  'Nnode=',this%Ncoord,&
  'nfault_elem=',this%nfault_elem

  !print*,'wave2fault=',this%wave2fault
  !print*,'fault2wave=',this%fault2wave
  !print*,shape(this%Tx0)
  !print*,'mpi_nn=',this%mpi_nn
  !print*,'mpi_ne=',this%mpi_ne
  !print*,'mpi_nmax=',this%mpi_nemax
  !print*,'coord=',this%coord
  !print*,'elem=',this%elem
  !print*,'elem=',this%part

end

subroutine write_mesh(this,filename,nproc)
  implicit none
  type(meshvar) :: this
  character(len=*) :: filename
  integer :: nproc

  !character(len=128) :: namestr
  !integer :: ios
  integer :: ierr,ncid,dimid(20),varid(50)
  !integer :: Nelem,Nnode

  !write(filename,'(a,i6.6,a)') 'data/meshVar',myrank,'.nc'
  ierr = nf90_create(trim(filename),NF90_WRITE,ncid)
  call check2(ierr,'nf90_create meshVar')
  !if (ierr /= nf90_noerr) stop 11
  !call check( nf90_create(trim(filename),NF90_WRITE,ncid) )

  ierr = nf90_def_dim(ncid,'Nnode',this%ncoord,dimid(5))
  call check2(ierr,'def_dim Nnode')
  ierr = nf90_def_dim(ncid,'Nelem',this%nelem,dimid(6))
  call check2(ierr,'def_dim Nelem')
  ierr = nf90_def_dim(ncid,'nfault_elem',this%nfault_elem,dimid(11))
  call check2(ierr,'def_dim nfault_elem')
  ierr = nf90_def_dim(ncid,'two',2,dimid(2))
  call check2(ierr,'def_dim two')
  ierr = nf90_def_dim(ncid,'three',3,dimid(3))
  call check2(ierr,'def_dim three')
  ierr = nf90_def_dim(ncid,'four',4,dimid(4))
  call check2(ierr,'def_dim four')

  if (nproc > 1) then
  ierr = nf90_def_dim(ncid,'mpi_nn',this%mpi_nn,dimid(7))
  call check2(ierr,'def_dim mpi_nn')
  ierr = nf90_def_dim(ncid,'mpi_ne',this%mpi_ne,dimid(8))
  call check2(ierr,'def_dim mpi_ne')
  ierr = nf90_def_dim(ncid,'mpi_nemax',this%mpi_nemax,dimid(9))
  call check2(ierr,'def_dim mpi_nemax')
  ierr = nf90_def_dim(ncid,'pinterfaces',this%pinterfaces,dimid(10))
  call check2(ierr,'def_dim pinterfaces')
  end if

  !call check( nf90_def_dim(ncid,'Nnode',this%ncoord,dimid(5)) )
  !call check( nf90_def_dim(ncid,'Nelem',this%nelem,dimid(6)) )
  !call check( nf90_def_dim(ncid,'nfault_elem',this%nfault_elem,dimid(11)) )
  !call check( nf90_def_dim(ncid,'mpi_nn',this%mpi_nn,dimid(7)) )
  !call check( nf90_def_dim(ncid,'mpi_ne',this%mpi_ne,dimid(8)) )
  !call check( nf90_def_dim(ncid,'mpi_nemax',this%mpi_nemax,dimid(9)) )
  !call check( nf90_def_dim(ncid,'pinterfaces',this%pinterfaces,dimid(10)) )
  !call check( nf90_def_dim(ncid,'two',2,dimid(2)) )
  !call check( nf90_def_dim(ncid,'three',3,dimid(3)) )
  !call check( nf90_def_dim(ncid,'four',4,dimid(4)) )

  ierr = nf90_def_var(ncid,'node',NF90_DOUBLE,  (/dimid(3),dimid(5)/),varid(1))
  call check2(ierr,'def_var node')
  ierr = nf90_def_var(ncid,'elem',NF90_INT,     (/dimid(4),dimid(6)/),varid(2))
  call check2(ierr,'def_var elem')
  ierr = nf90_def_var(ncid,'neighbor',NF90_INT, (/dimid(4),dimid(6)/),varid(3))
  call check2(ierr,'def_var neighbor')
  ierr = nf90_def_var(ncid,'face',NF90_INT,     (/dimid(4),dimid(6)/),varid(4))
  call check2(ierr,'def_var face')
  ierr = nf90_def_var(ncid,'direction',NF90_INT,(/dimid(4),dimid(6)/),varid(5))
  call check2(ierr,'def_var direction')
  ierr = nf90_def_var(ncid,'bctype',NF90_INT,   (/dimid(4),dimid(6)/),varid(6))
  call check2(ierr,'def_var bctype')
  ierr = nf90_def_var(ncid,'fluxtype',NF90_INT,   (/dimid(4),dimid(6)/),varid(30))
  call check2(ierr,'def_var fluxtype')
  ierr = nf90_def_var(ncid,'elemtype',NF90_INT,   dimid(6),varid(7))
  call check2(ierr,'def_var elemtype')
  ierr = nf90_def_var(ncid,'rho',NF90_DOUBLE,dimid(6),varid(8))
  call check2(ierr,'def_var rho')
  ierr = nf90_def_var(ncid,'vp' ,NF90_DOUBLE,dimid(6),varid(9))
  call check2(ierr,'def_var vp')
  ierr = nf90_def_var(ncid,'vs' ,NF90_DOUBLE,dimid(6),varid(10))
  call check2(ierr,'def_var vs')

  ierr = nf90_def_var(ncid,'Tx0',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(20))
  call check2(ierr,'def_var Tx0')
  ierr = nf90_def_var(ncid,'Ty0',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(21))
  call check2(ierr,'def_var Ty0')
  ierr = nf90_def_var(ncid,'Tz0',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(22))
  call check2(ierr,'def_var Tz0')
  ierr = nf90_def_var(ncid,'dTx0',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(33))
  call check2(ierr,'def_var dTx0')
  ierr = nf90_def_var(ncid,'dTy0',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(34))
  call check2(ierr,'def_var dTy0')
  ierr = nf90_def_var(ncid,'dTz0',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(35))
  call check2(ierr,'def_var dTz0')
  ierr = nf90_def_var(ncid,'mu_s',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(23))
  call check2(ierr,'def_var mu_s')
  ierr = nf90_def_var(ncid,'mu_d',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(24))
  call check2(ierr,'def_var mu_d')
  ierr = nf90_def_var(ncid,'Dc',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(25))
  call check2(ierr,'def_var Dc')
  ierr = nf90_def_var(ncid,'C0',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(26))
  call check2(ierr,'def_var C0')
  ! rate state
  ierr = nf90_def_var(ncid,'a',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(27))
  call check2(ierr,'def_var a')
  ierr = nf90_def_var(ncid,'b',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(28))
  call check2(ierr,'def_var b')
  ierr = nf90_def_var(ncid,'Vw',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(31))
  call check2(ierr,'def_var Vw')
  ierr = nf90_def_var(ncid,'state',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(29))
  call check2(ierr,'def_var state')
  ! thermal pressurization
  ierr = nf90_def_var(ncid,'TP_hy',NF90_DOUBLE,(/dimid(3),dimid(4),dimid(11)/),varid(32))
  call check2(ierr,'def_var TP_hy')

  if (nproc > 1) then
  ierr = nf90_def_var(ncid,'mpi_rho',NF90_DOUBLE,(/dimid(10),dimid(3)/),varid(15))
  call check2(ierr,'def_var mpi_rho')
  ierr = nf90_def_var(ncid,'mpi_vp',NF90_DOUBLE,(/dimid(10),dimid(3)/),varid(16))
  call check2(ierr,'def_var mpi_vp')
  ierr = nf90_def_var(ncid,'mpi_vs',NF90_DOUBLE,(/dimid(10),dimid(3)/),varid(17))
  call check2(ierr,'def_var mpi_vs')

  ierr = nf90_def_var(ncid,'mpi_neighbor',NF90_INT,dimid(7),varid(11))
  call check2(ierr,'def_var mpi_neighbor')
  ierr = nf90_def_var(ncid,'mpi_connection',NF90_INT,(/dimid(7),dimid(8),dimid(2)/),varid(12))
  call check2(ierr,'def_var mpi_connnetion')
  ierr = nf90_def_var(ncid,'mpi_ibool',NF90_INT,(/dimid(4),dimid(9)/),varid(13))
  call check2(ierr,'def_var mpi_ibool')
  ierr = nf90_def_var(ncid,'mpi_interface',NF90_INT,(/dimid(4),dimid(4),dimid(6)/),varid(14))
  call check2(ierr,'def_var mpi_interface')
  end if

  ierr = nf90_enddef(ncid)

  ierr = nf90_put_var(ncid,varid(1),this%coord)
  call check2(ierr,'put_var coord')
  ierr = nf90_put_var(ncid,varid(2),this%elem)
  call check2(ierr,'put_var elem')
  ierr = nf90_put_var(ncid,varid(3),this%neighbor)
  call check2(ierr,'put_var neighbor')
  ierr = nf90_put_var(ncid,varid(4),this%face)
  call check2(ierr,'put_var face')
  ierr = nf90_put_var(ncid,varid(5),this%direction)
  call check2(ierr,'put_var direction')
  ierr = nf90_put_var(ncid,varid(6),this%bctype)
  call check2(ierr,'put_var bctype')
  ierr = nf90_put_var(ncid,varid(30),this%fluxtype)
  call check2(ierr,'put_var fluxtype')
  ierr = nf90_put_var(ncid,varid(7),this%elemtype)
  call check2(ierr,'put_var elemtype')
  ierr = nf90_put_var(ncid,varid(8),this%rho)
  call check2(ierr,'put_var rho')
  ierr = nf90_put_var(ncid,varid(9),this%vp)
  call check2(ierr,'put_var vp')
  ierr = nf90_put_var(ncid,varid(10),this%vs)
  call check2(ierr,'put_var vs')
  ierr = nf90_put_var(ncid,varid(15),this%mpi_rho)
  if (nproc > 1) then
  call check2(ierr,'put_var mpi_rho')
  ierr = nf90_put_var(ncid,varid(16),this%mpi_vp)
  call check2(ierr,'put_var mpi_vp')
  ierr = nf90_put_var(ncid,varid(17),this%mpi_vs)
  call check2(ierr,'put_var mpi_vs')
  ierr = nf90_put_var(ncid,varid(11),this%mpi_neighbor)
  call check2(ierr,'put_var mpi_neighbor')
  ierr = nf90_put_var(ncid,varid(12),this%mpi_connection)
  call check2(ierr,'put_var mpi_connection')
  ierr = nf90_put_var(ncid,varid(13),this%mpi_ibool)
  call check2(ierr,'put_var ibool')
  ierr = nf90_put_var(ncid,varid(14),this%mpi_interface)
  call check2(ierr,'put_var mpi_interface')
  end if

  ierr = nf90_put_var(ncid,varid(20),this%Tx0)
  call check2(ierr,'put_var Tx0')
  ierr = nf90_put_var(ncid,varid(21),this%Ty0)
  call check2(ierr,'put_var Ty0')
  ierr = nf90_put_var(ncid,varid(22),this%Tz0)
  call check2(ierr,'put_var Tz0')
  ierr = nf90_put_var(ncid,varid(33),this%dTx0)
  call check2(ierr,'put_var dTx0')
  ierr = nf90_put_var(ncid,varid(34),this%dTy0)
  call check2(ierr,'put_var dTy0')
  ierr = nf90_put_var(ncid,varid(35),this%dTz0)
  call check2(ierr,'put_var dTz0')
  ierr = nf90_put_var(ncid,varid(23),this%mu_s)
  call check2(ierr,'put_var mu_s')
  ierr = nf90_put_var(ncid,varid(24),this%mu_d)
  call check2(ierr,'put_var mu_d')
  ierr = nf90_put_var(ncid,varid(25),this%Dc)
  call check2(ierr,'put_var Dc')
  ierr = nf90_put_var(ncid,varid(26),this%C0)
  call check2(ierr,'put_var C0')
  ! rate state
  ierr = nf90_put_var(ncid,varid(27),this%a)
  call check2(ierr,'put_var a')
  ierr = nf90_put_var(ncid,varid(28),this%b)
  call check2(ierr,'put_var b')
  ierr = nf90_put_var(ncid,varid(31),this%Vw)
  call check2(ierr,'put_var Vw')
  ierr = nf90_put_var(ncid,varid(29),this%state)
  call check2(ierr,'put_var state')
  ! thermal pressurization
  ierr = nf90_put_var(ncid,varid(32),this%TP_hy)
  call check2(ierr,'put_var TP_hy')

  ierr = nf90_close(ncid)
  call check2(ierr,'nf90_close meshVar')
end subroutine

end program
