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

module mod_types

  use mod_para, only : RKIND

  implicit none

  public :: meshvar

  type :: meshvar
    integer :: nelem,ncoord
    integer :: nglob

    ! mpi
    integer :: rank,nproc
    integer :: mpi_nn
    integer :: mpi_ne
    integer :: mpi_nemax
    integer :: pinterfaces

    integer :: nrecv
    integer :: body_nrecv
    integer :: nfault_elem
    integer :: nfault_face ! one element may have >1 fault face
    integer :: nfree_face
    integer :: irk,nrk

    real(kind=RKIND) :: dtfactor,deltat
    real(kind=RKIND) :: current_time
    real(kind=RKIND) :: rmin,rmax
    real(kind=RKIND) :: xmin,xmax
    real(kind=RKIND) :: ymin,ymax
    real(kind=RKIND) :: zmin,zmax

    real(kind=RKIND),dimension(:,:),pointer :: coord=>null()
    integer,dimension(:,:),pointer :: elem => null()
    integer,dimension(:,:),pointer :: neigh => null()
    integer,dimension(:,:),pointer :: face => null()
    integer,dimension(:,:),pointer :: direction => null()
    integer,dimension(:,:),pointer :: bctype => null()
    integer,dimension(:,:),pointer :: fluxtype => null()
    integer,dimension(:),pointer :: elemtype => null()
    !real(kind=RKIND),dimension(:,:),allocatable :: coord
    !integer,dimension(:,:),allocatable :: elem
    !integer,dimension(:,:),allocatable :: neigh
    !integer,dimension(:,:),allocatable :: face
    !integer,dimension(:,:),allocatable :: direction
    !integer,dimension(:,:),allocatable :: bctype
    !integer,dimension(:),allocatable :: elemtype
    !integer :: vmapM(Nfp,Nfaces)
    integer,dimension(:,:,:),pointer :: vmapM => null()
    integer,dimension(:,:,:),pointer :: vmapP => null()
    !integer :: flipped_index(Nfp,6)
    integer,dimension(:,:),pointer :: flipped_index => null()

    ! mpi
    integer,dimension(:),pointer :: mpi_neighbor => null()
    integer,dimension(:,:,:),pointer :: mpi_connection => null()
    integer,dimension(:,:),pointer :: mpi_ibool => null()
    integer,dimension(:,:,:),pointer :: mpi_interface => null()
    real(kind=rkind),dimension(:,:),pointer :: mpi_vp => null()
    real(kind=rkind),dimension(:,:),pointer :: mpi_vs => null()
    real(kind=rkind),dimension(:,:),pointer :: mpi_rho => null()
    !integer,dimension(:),allocatable :: mpi_neighbor
    !integer,dimension(:,:,:),allocatable :: mpi_connection
    !integer,dimension(:,:),allocatable :: mpi_ibool
    !integer,dimension(:,:,:),allocatable :: mpi_interface

    ! fault
    integer,dimension(:),pointer :: fault2wave => null()
    integer,dimension(:),pointer :: wave2fault => null()
    !
    real(kind=rkind),dimension(:,:,:),pointer :: tau0n => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tau0m => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tau0l => null()
    real(kind=rkind),dimension(:,:,:),pointer :: dtau0n => null()
    real(kind=rkind),dimension(:,:,:),pointer :: dtau0m => null()
    real(kind=rkind),dimension(:,:,:),pointer :: dtau0l => null()
    real(kind=rkind),dimension(:,:,:),pointer :: stress => null()
    real(kind=rkind),dimension(:,:,:),pointer :: stress1 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: stress2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sigma => null()
    real(kind=rkind),dimension(:,:,:),pointer :: slip => null()
    real(kind=rkind),dimension(:,:,:),pointer :: slip1 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: slip2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mslip1 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tslip1 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mslip2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tslip2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sliprate => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sliprate1 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sliprate2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mu_s => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mu_d => null()
    real(kind=rkind),dimension(:,:,:),pointer :: Dc => null()
    real(kind=rkind),dimension(:,:,:),pointer :: C0 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: ruptime => null()

    ! rate state
    real(kind=rkind),dimension(:,:,:),pointer :: a => null()
    real(kind=rkind),dimension(:,:,:),pointer :: b => null()
    real(kind=rkind),dimension(:,:,:),pointer :: Vw => null()
    real(kind=rkind),dimension(:,:,:),pointer :: state => null()
    real(kind=rkind),dimension(:,:,:),pointer :: hstate => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mstate => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tstate => null()
    ! thermal pressurization
    integer :: TP_n
    real(kind=rkind),dimension(:,:,:,:),pointer :: TP_Theta => null()
    real(kind=rkind),dimension(:,:,:,:),pointer :: TP_Sigma => null()
    real(kind=rkind),dimension(:,:,:),pointer :: TP_T => null()
    real(kind=rkind),dimension(:,:,:),pointer :: TP_P => null()
    real(kind=rkind),dimension(:,:,:),pointer :: TP_hy => null()
    real(kind=rkind),dimension(:),pointer :: DFinv => null()
    real(kind=rkind),dimension(:),pointer :: Dwn => null()


    ! matrices
    !real(kind=rkind),dimension(Np,Np) :: Vdm,invVdm,Mass,Dr,Ds,Dt!,Filter
    !integer :: Fmask(Nfp,Nfaces)
    !real(kind=rkind) :: LIFT(Np,Nfaces*Nfp)
    integer,dimension(:,:),pointer :: Fmask => null()
    real(kind=rkind),dimension(:,:),pointer :: Vdm => null()
    real(kind=rkind),dimension(:,:),pointer :: invVdm => null()
    real(kind=rkind),dimension(:,:),pointer :: Mass => null()
    real(kind=rkind),dimension(:,:),pointer :: Dr => null()
    real(kind=rkind),dimension(:,:),pointer :: Ds => null()
    real(kind=rkind),dimension(:,:),pointer :: Dt => null()
    real(kind=rkind),dimension(:,:),pointer :: LIFT => null()


    real(kind=rkind),dimension(:),pointer :: r => null()
    real(kind=rkind),dimension(:),pointer :: s => null()
    real(kind=rkind),dimension(:),pointer :: t => null()
    real(kind=rkind),dimension(:),pointer :: vx => null()
    real(kind=rkind),dimension(:),pointer :: vy => null()
    real(kind=rkind),dimension(:),pointer :: vz => null()
    real(kind=rkind),dimension(:),pointer :: rx => null()
    real(kind=rkind),dimension(:),pointer :: ry => null()
    real(kind=rkind),dimension(:),pointer :: rz => null()
    real(kind=rkind),dimension(:),pointer :: sx => null()
    real(kind=rkind),dimension(:),pointer :: sy => null()
    real(kind=rkind),dimension(:),pointer :: sz => null()
    real(kind=rkind),dimension(:),pointer :: tx => null()
    real(kind=rkind),dimension(:),pointer :: ty => null()
    real(kind=rkind),dimension(:),pointer :: tz => null()
    real(kind=rkind),dimension(:),pointer :: jac => null()
    real(kind=rkind),dimension(:,:),pointer :: nx => null()
    real(kind=rkind),dimension(:,:),pointer :: ny => null()
    real(kind=rkind),dimension(:,:),pointer :: nz => null()
    real(kind=rkind),dimension(:,:),pointer :: mx => null()
    real(kind=rkind),dimension(:,:),pointer :: my => null()
    real(kind=rkind),dimension(:,:),pointer :: mz => null()
    real(kind=rkind),dimension(:,:),pointer :: lx => null()
    real(kind=rkind),dimension(:,:),pointer :: ly => null()
    real(kind=rkind),dimension(:,:),pointer :: lz => null()
    real(kind=rkind),dimension(:,:),pointer :: sJ => null()
    real(kind=rkind),dimension(:,:),pointer :: Fscale => null()

    ! media
    real(kind=rkind),dimension(:),pointer :: vp => null()
    real(kind=rkind),dimension(:),pointer :: vs => null()
    real(kind=rkind),dimension(:),pointer :: rho => null()
    real(kind=rkind),dimension(:),pointer :: zp => null()
    real(kind=rkind),dimension(:),pointer :: zs => null()
    real(kind=rkind),dimension(:),pointer :: lam => null()
    real(kind=rkind),dimension(:),pointer :: mu => null()
    !real(kind=rkind),dimension(:),allocatable :: vp
    !real(kind=rkind),dimension(:),allocatable :: vs
    !real(kind=rkind),dimension(:),allocatable :: rho
    !real(kind=rkind),dimension(:),allocatable :: zp
    !real(kind=rkind),dimension(:),allocatable :: zs
    !real(kind=rkind),dimension(:),allocatable :: lam
    !real(kind=rkind),dimension(:),allocatable :: mu
    real(kind=rkind),dimension(:),pointer :: damp => null()
  !@
  !@  ! fault recvs
  !@  integer,dimension(:),pointer :: recv_fid => null()
  !@  integer,dimension(:),pointer :: recv_i => null()
  !@  integer,dimension(:),pointer :: recv_ie => null()
  !@  real(kind=rkind),dimension(:),pointer :: recv_refx => null()
  !@  real(kind=rkind),dimension(:,:,:),pointer :: recv_buffer => null()
  !@  real(kind=rkind),dimension(:,:,:),pointer :: fault_buffer => null()
  !@
  !@  real(kind=rkind),dimension(:,:,:),pointer :: surface_buffer => null()
  !@
  !@  ! body recvs
  !@  integer,dimension(:),pointer :: body_recv_fid => null()
  !@  integer,dimension(:),pointer :: body_recv_i => null()
  !@  integer,dimension(:),pointer :: body_recv_j => null()
  !@  integer,dimension(:),pointer :: body_recv_ie => null()
  !@  real(kind=rkind),dimension(:),pointer :: body_recv_refx => null()
  !@  real(kind=rkind),dimension(:),pointer :: body_recv_refy => null()
    ! recvs
    integer,dimension(:),pointer :: recv_id => null()
    integer,dimension(:),pointer :: recv_bctype => null()
    integer,dimension(:),pointer :: recv_elem => null()
    integer,dimension(:),pointer :: recv_face => null()
    real(kind=rkind),dimension(:,:),pointer :: recv_coord => null()
    real(kind=rkind),dimension(:,:),pointer :: recv_normal => null()
    real(kind=rkind),dimension(:,:,:),pointer :: recv_buffer => null()

  end type

  type wavefield
    ! u(Np*Nelem,9)
    real(kind=RKIND),allocatable,dimension(:,:) :: u,hu,tu,mu
    real(kind=RKIND),allocatable,dimension(:,:) :: displ
  end type

  type fault
    ! slip(Nfp,Nface,Nelem_fault)
    !real(kind=RKIND),allocatable,dimension(:,:) :: &
    !slip,&
    !hslip,&
    !mslip,&
    !tslip
    real(kind=rkind),dimension(:,:,:),pointer :: tau0n => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tau0m => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tau0l => null()
    real(kind=rkind),dimension(:,:,:),pointer :: dtau0n => null()
    real(kind=rkind),dimension(:,:,:),pointer :: dtau0m => null()
    real(kind=rkind),dimension(:,:,:),pointer :: dtau0l => null()
    real(kind=rkind),dimension(:,:,:),pointer :: stress => null()
    real(kind=rkind),dimension(:,:,:),pointer :: stress1 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: stress2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sigma => null()
    real(kind=rkind),dimension(:,:,:),pointer :: slip => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mslip => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tslip => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sliprate => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sliprate1 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sliprate2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mu_s => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mu_d => null()
    real(kind=rkind),dimension(:,:,:),pointer :: Dc => null()
    real(kind=rkind),dimension(:,:,:),pointer :: C0 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: ruptime => null()

    ! rate state
    real(kind=rkind),dimension(:,:,:),pointer :: a => null()
    real(kind=rkind),dimension(:,:,:),pointer :: b => null()
    real(kind=rkind),dimension(:,:,:),pointer :: Vw => null()
    real(kind=rkind),dimension(:,:,:),pointer :: state => null()
    real(kind=rkind),dimension(:,:,:),pointer :: hstate => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mstate => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tstate => null()
    ! thermal pressurization
    integer :: TP_n
    real(kind=rkind),dimension(:,:,:,:),pointer :: TP_Theta => null()
    real(kind=rkind),dimension(:,:,:,:),pointer :: TP_Sigma => null()
    real(kind=rkind),dimension(:,:,:),pointer :: TP_T => null()
    real(kind=rkind),dimension(:,:,:),pointer :: TP_P => null()
    real(kind=rkind),dimension(:,:,:),pointer :: TP_hy => null()
    real(kind=rkind),dimension(:),pointer :: DFinv => null()
    real(kind=rkind),dimension(:),pointer :: Dwn => null()

  end type

  type buffvar
    real(kind=rkind), dimension(:,:), allocatable :: q_send,q_rec
    real(kind=rkind), dimension(:,:,:,:), allocatable :: qi
  end type

  type mytimer
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    integer,dimension(8) :: values
  end type

contains

end module
