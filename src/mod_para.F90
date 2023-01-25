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

module mod_para

! default pOrder
#ifndef pOrder
#define pOrder 2
#endif

  implicit none

  character(len=256) :: problem
  character(len=256) :: mesh_dir
  character(len=256) :: data_dir

  integer,parameter :: SIZE_REAL = 4
  integer,parameter :: SIZE_DOUBLE = 8
  integer,parameter :: RKIND = SIZE_DOUBLE
  !integer,parameter :: CUSTOM_REAL = SIZE_DOUBLE
  integer,parameter :: CUSTOM_REAL = RKIND
  integer,parameter :: MAX_NUM_RECV = 1000

  !integer,parameter :: nsurface = 4
  integer,parameter :: Nfaces = 4
  integer,parameter :: Order = pOrder
  integer,parameter :: NGLL = Order+1
  integer,parameter :: Nfp = (order+1)*(order+2)/2
  integer,parameter :: Np = (Order+1)*(Order+2)*(Order+3)/6
  integer,parameter :: Nvar = 9
  integer,parameter :: dimens = 9

  real(kind=RKIND),parameter :: PI = 3.141592653589793238463d0
  real(kind=RKIND),parameter :: SQRT2PI = 2.506628274631000d0
  real(kind=RKIND),parameter :: zero = 0.0
  real(kind=RKIND),parameter :: one = 1.0
  real(kind=RKIND),parameter :: two = 2.0
  real(kind=RKIND),parameter :: EPS = 1e-5

  integer,parameter :: BC_IN    = 0
  integer,parameter :: BC_FREE  = 1
  integer,parameter :: BC_OUT   = 2
  integer,parameter :: BC_FAULT = 100

  ! 0: classic RK4
  ! 1: RK54
  ! 2: RK33
  ! 3: RK22
  integer,parameter :: timeIntegrationMethod = 1

  ! Low storage Runge-Kutta coefficients

  ! Heun (2,2)
  real(kind=RKIND),parameter :: rk2a(2) = (/ 0d0, -1d0 /)
  real(kind=RKIND),parameter :: rk2b(2) = (/ 1d0, 0.5d0 /)
  real(kind=RKIND),parameter :: rk2c(0:2) = (/0d0, 1d0, 1d0 /)

  ! Williamson (3,3)
  real(kind=RKIND),parameter :: rk3a(3) = (/ 0d0, -5d0/9d0, -153d0/128d0 /)
  real(kind=RKIND),parameter :: rk3b(3) = (/ 1d0/3d0, 15d0/16d0, 8d0/15d0 /)
  real(kind=RKIND),parameter :: rk3c(0:3) = (/0d0, 1d0/3d0, 3d0/4d0, 1d0 /)

  ! Kennedy-Carpenter (5,4)
  real(kind=RKIND),parameter :: rk4a(5) = (/ &
          0.0, &
          -567301805773.0/1357537059087.0, &
          -2404267990393.0/2016746695238.0, &
          -3550918686646.0/2091501179385.0, &
          -1275806237668.0/842570457699.0/)
  real(kind=RKIND),parameter :: rk4b(5) = (/ &
          1432997174477.0/9575080441755.0, &
          5161836677717.0/13612068292357.0, &
          1720146321549.0/2090206949498.0, &
          3134564353537.0/4481467310338.0, &
          2277821191437.0/14882151754819.0/)
  real(kind=RKIND),parameter :: rk4c(0:5) = (/ &
          0.0, &
          1432997174477.0/9575080441755.0, &
          2526269341429.0/6820363962896.0, &
          2006345519317.0/3224310063776.0, &
          2802321613138.0/2924317926251.0, &
          1.0 /)

  ! mpi
  logical :: masternode

  logical :: nice_print

  real(kind=RKIND) :: simu_time_max
  real(kind=RKIND) :: cfl_number
  real(kind=RKIND) :: timestep

  integer :: fault_snap_skip
  integer :: grdsurf_snap_skip

  integer :: flux_method

  integer :: use_damp
  integer :: use_pml

  ! Friction laws
  ! 0 : linear slip weakening (default)
  ! 1 : rate state, ageing law
  ! 2 : rate state, slip law
  ! 3 : rate state, slip law, flash heating
  ! 4 : time weakening
  integer :: friction_law

  real(kind=rkind) :: RS_f0
  real(kind=rkind) :: RS_V0
  real(kind=rkind) :: RS_fw

  ! off-fault plasticity
  integer :: plasticity
  real(kind=rkind) :: cohesion
  real(kind=rkind) :: blkfric
  real(kind=rkind) :: Tvisc
  real(kind=rkind) :: coef_byy
  real(kind=rkind) :: coef_bxx
  real(kind=rkind) :: coef_bxy
  real(kind=rkind) :: fluidpres_profile_h1
  real(kind=rkind) :: fluidpres_profile_h2
  real(kind=rkind) :: fluidpres_profile_o1
  real(kind=rkind) :: fluidpres_profile_o2

  ! thermal pressurization
  integer :: thermalpressure

  ! io
  integer :: export_grdsurf_velo
  integer :: export_grdsurf_displ
  integer :: export_grdsurf_strain
  integer :: export_media
  integer :: export_wave
  integer :: export_wave_component
  integer :: export_wave_timestep

  ! smoothly loading stress perturbation: T0+coef*dT0
  integer :: smooth_load
  real(kind=RKIND) :: smooth_load_time

  ! Initial condition
  integer :: initial_condition_wave
  real(kind=rkind) :: src_loc(3)
  real(kind=rkind) :: src_gaussian_width
  real(kind=rkind) :: src_mxx
  real(kind=rkind) :: src_myy
  real(kind=rkind) :: src_mzz
  real(kind=rkind) :: src_myz
  real(kind=rkind) :: src_mxz
  real(kind=rkind) :: src_mxy
  real(kind=rkind) :: src_m0

  ! parameters for time weakening law
  real(kind=rkind) :: nucleate_y0
  real(kind=rkind) :: nucleate_z0
  real(kind=rkind) :: nucleate_rcrit
  real(kind=rkind) :: nucleate_Vrup
  real(kind=rkind) :: TimeForcedRup


  ! these variables will be used in each time step:
  ! variables in each element
  !real(kind=RKIND),dimension(:,:),allocatable :: &
  !  Ue,F1,F2,F3,fluxs
  !real(kind=RKIND),dimension(:),allocatable :: &
  !  dF1dr,dF1ds,dF1dt, &
  !  dF2dr,dF2ds,dF2dt, &
  !  dF3dr,dF3ds,dF3dt, &
  !  flux
  !real(kind=RKIND),dimension(Np,Nvar) :: Ue,F1,F2,F3 ! (Np,Nvar)
  !real(kind=RKIND),dimension(Np) :: dF1dr,dF1ds,dF1dt
  !real(kind=RKIND),dimension(Np) :: dF2dr,dF2ds,dF2dt
  !real(kind=RKIND),dimension(Np) :: dF3dr,dF3ds,dF3dt
  !real(kind=RKIND),dimension(Nfp*Nfaces,Nvar) :: fluxs 
  !real(kind=RKIND),dimension(Nvar) :: flux 
  !real(kind=RKIND),dimension(:,:),pointer :: fluxes=>null()

!contains


end module
