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

PROGRAM MAIN

use mod_para,       only : CUSTOM_REAL,           &
                           RKIND,                 &
                           Nvar, Np, Nfaces, Nfp, &
                           masternode,            &
                           use_damp,              &
                           timestep,              &
                           simu_time_max,         &
                           cfl_number,            &
                           fault_snap_skip,       &
                           grdsurf_snap_skip,     &
                           thermalpressure,       &
                           plasticity,            &
                           nice_print,            &
                           rk2a, rk2b,            &
                           rk3a, rk3b,            &
                           rk4a, rk4b,            &
                           timeIntegrationMethod, &
                           export_wave,           &
                           export_wave_timestep
use mod_read,       only : read_parameters
use mod_types,      only : meshvar,               &
                           buffvar,               &
                           mytimer
use mod_mesh,       only : readMeshVar
use mod_geometry,   only : build_geometry
use mod_wave,       only : rhs,                   &
                           update_info,           &
                           init_wave,             &
                           write_wave_vtk
use mod_plastic,    only : update_plastic
use mod_thermpress, only : init_thermpress,       &
                           add_TP
use mod_fault,      only : fault_init
use mod_io_fault,   only : fault_io_init,         &
                           fault_io_save,         &
                           fault_io_end
use mod_io_grdsurf, only : grdsurf_io_init,       &
                           grdsurf_io_save,       &
                           grdsurf_io_end,        &
                           update_pgv
use mod_io_recv,    only : recv_io_init,          &
                           recv_io_save,          &
                           recv_io_end
use mod_io_energy,  only : energy_io_init,          &
                           energy_io_save,          &
                           energy_io_end
!use mod_source
use mod_recv,       only : locate_recvs

use mpi
use mod_mpi,        only : sync_print,            &
                           minval_real_all
use mod_exchange,   only : init_buff,             &
                           exchange_data
use mod_damp,       only : init_damp,             &
                           cal_coord_range
use mod_check,      only : check_geometry

implicit none

type(meshvar) :: mesh
type(buffvar) :: buff
real(kind=RKIND),allocatable,dimension(:,:) :: u,hu,tu,mu
real(kind=RKIND),allocatable,dimension(:,:) :: displ
integer :: myrank,nproc,ierr
integer :: i,it,nt,irk
real(kind=RKIND) :: dt,dtfactor
double precision :: wtime1,wtime2
real :: time_start,time_end
type(mytimer) :: timer1

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
mesh%rank = myrank
mesh%nproc = nproc

if (myrank==0) then
  masternode = .true.
else
  masternode = .false.
end if

nice_print = .true.

call read_parameters()

!call get_command_argument(1,data_dir)
!if (len_trim(data_dir)==0) then
!  data_dir = 'data'
!end if
!if(myrank==0) &
!write(*,'(a,a)') 'data_dir = ',trim(data_dir)

call sync_print()

call readMeshVar(mesh)

call sync_print()

print*, 'rank=',myrank,'ncoord=',mesh%ncoord,'nelem=',mesh%nelem

call sync_print()

write(*,'(a,i8,a,2f10.4,a,2f10.4,a,2f10.4)') &
    ' rank=',myrank,&
    ' vp=',minval(mesh%vp),maxval(mesh%vp),&
    ' vs=',minval(mesh%vs),maxval(mesh%vs),&
    ' rho=',minval(mesh%rho),maxval(mesh%rho)

call sync_print()

call build_geometry(mesh)

! check
!call check_geometry(mesh)

call sync_print()

write(*,'(a,i8,a,2f10.4,a,f10.4)') &
    ' rank=',myrank,&
    ' radius=',mesh%rmin,mesh%rmax,&
    ' dtfactor=',mesh%dtfactor

call sync_print()

call locate_recvs(mesh)

allocate(u(Np*mesh%Nelem,Nvar))
allocate(hu(Np*mesh%Nelem,Nvar))
allocate(tu(Np*mesh%Nelem,Nvar))
allocate(mu(Np*mesh%Nelem,Nvar))
allocate(displ(Np*mesh%Nelem,3))
tu=0;mu=0;hu=0;u=0

!allocate(fluxes(mesh%nelem*Nfaces*Nfp,Nvar))

call init_wave(mesh,u)

call cal_coord_range(mesh)

call init_damp(mesh)

call fault_init(mesh)

call init_thermpress(mesh)

print*,'rank=',myrank,&
'nfault_elem=',mesh%nfault_elem,&
'nfault_face=',mesh%nfault_face,&
'nfree_face=',mesh%nfree_face

call sync_print()

call fault_io_init(mesh)

call grdsurf_io_init(mesh)

call recv_io_init(mesh)

call energy_io_init(mesh)

call sync_print()


! init mpi buffer
call init_buff(mesh,buff)

call minval_real_all(mesh%dtfactor,dtfactor,CUSTOM_REAL)

if(myrank==0) print*,'reduced dtfactor = ',dtfactor

call sync_print()

dt = cfl_number*dtfactor

mesh%deltat = dt

if (timestep > 0) then
  if(masternode) print*,'timestep=',mesh%deltat,'now set to',timestep
  dt = timestep
  mesh%deltat = dt
end if

nt = int(simu_time_max/dt)+1

if(myrank==0) print*,'dt=',dt,'nt=',nt

! RK method
if(timeIntegrationMethod == 0) then
  mesh%nrk = 4
elseif(timeIntegrationMethod == 1) then
  mesh%nrk = 5
elseif(timeIntegrationMethod == 2) then
  mesh%nrk = 3
elseif(timeIntegrationMethod == 3) then
  mesh%nrk = 2
end if

displ = 0

call cpu_time(time_start)
wtime1 = MPI_WTIME()
do it = 1,nt

  mesh%current_time = (it-1)*dt

  mu = u
  mesh%mslip1 = mesh%slip1
  mesh%mslip2 = mesh%slip2
  tu = 0.0
  mesh%tslip1 = 0.0
  mesh%tslip2 = 0.0
  ! rate state
  mesh%mstate = mesh%state
  mesh%tstate = 0.0

  if (mod(it-1,fault_snap_skip) == 0) then
    call fault_io_save(mesh,(it-1)/fault_snap_skip+1)
    wtime2 = MPI_WTIME()
    if(myrank==0) then
      call date_and_time(timer1%date,timer1%time,timer1%zone,timer1%values)
      write(*,'(a,1x,a,1x,a,1x,a,i8,1x,f10.4,a)') &
      timer1%date,timer1%time,timer1%zone, &
      '@ timestep ',it,(wtime2-wtime1)/dble(fault_snap_skip),' s/step'
    end if
    wtime1 = wtime2
  end if

  call update_pgv(mesh,u)
  if (mod(it-1,grdsurf_snap_skip) == 0) then
    call grdsurf_io_save(mesh,u,displ,(it-1)/grdsurf_snap_skip+1)
  end if

  call recv_io_save(mesh,u,it)

  do irk = 1,mesh%nrk
    !tu=0
    mesh%irk = irk

    call exchange_data(mesh,u,buff)

    ! smooth
    !mu = u
    !call avg_face(mesh,mu,u,qi)

    call rhs(mesh,u,buff%qi,hu)

    if (timeIntegrationMethod == 0) then
      if(irk==1)then
        u  = mu+0.5d0*dt*hu
        tu = mu+1.0d0/6.0d0*dt*hu

        mesh%slip1  = mesh%mslip1+0.5d0*dt*mesh%sliprate1
        mesh%tslip1 = mesh%mslip1+1.0d0/6.0d0*dt*mesh%sliprate1
        mesh%slip2  = mesh%mslip2+0.5d0*dt*mesh%sliprate2
        mesh%tslip2 = mesh%mslip2+1.0d0/6.0d0*dt*mesh%sliprate2
        ! rate state
        mesh%state  = mesh%mstate+0.5d0*dt*mesh%hstate
        mesh%tstate = mesh%mstate+1.0d0/6.0d0*dt*mesh%hstate
      elseif(irk==2)then
        u  = mu+0.5d0*dt*hu
        tu = tu+1.0d0/3.0d0*dt*hu

        mesh%slip1  = mesh%mslip1+0.5d0*dt*mesh%sliprate1
        mesh%tslip1 = mesh%tslip1+1.0d0/3.0d0*dt*mesh%sliprate1
        mesh%slip2  = mesh%mslip2+0.5d0*dt*mesh%sliprate2
        mesh%tslip2 = mesh%tslip2+1.0d0/3.0d0*dt*mesh%sliprate2
        ! rate state
        mesh%state  = mesh%mstate+0.5d0*dt*mesh%hstate
        mesh%tstate = mesh%tstate+1.0d0/3.0d0*dt*mesh%hstate
      elseif(irk==3)then
        u  = mu+1.0d0*dt*hu
        tu = tu+1.0d0/3.0d0*dt*hu

        mesh%slip1  = mesh%mslip1+1.0d0*dt*mesh%sliprate1
        mesh%tslip1 = mesh%tslip1+1.0d0/3.0d0*dt*mesh%sliprate1
        mesh%slip2  = mesh%mslip2+1.0d0*dt*mesh%sliprate2
        mesh%tslip2 = mesh%tslip2+1.0d0/3.0d0*dt*mesh%sliprate2
        ! rate state
        mesh%state  = mesh%mstate+1.0d0*dt*mesh%hstate
        mesh%tstate = mesh%tstate+1.0d0/3.0d0*dt*mesh%hstate
      elseif(irk==4)then
        u  = tu+1.0d0/6.0d0*dt*hu

        mesh%slip1 = mesh%tslip1+1.0d0/6.0d0*dt*mesh%sliprate1
        mesh%slip2 = mesh%tslip2+1.0d0/6.0d0*dt*mesh%sliprate2
        ! rate state
        mesh%state  = mesh%tstate+1.0d0/6.0d0*dt*mesh%hstate
      endif

    else if (timeIntegrationMethod == 1) then

      tu = rk4a(irk)*tu + dt*hu
      u = u + rk4b(irk)*tu

      mesh%tslip1 = rk4a(irk)*mesh%tslip1 + dt*mesh%sliprate1
      mesh%slip1 = mesh%slip1 + rk4b(irk)*mesh%tslip1
      mesh%tslip2 = rk4a(irk)*mesh%tslip2 + dt*mesh%sliprate2
      mesh%slip2 = mesh%slip2 + rk4b(irk)*mesh%tslip2
      ! rate state
      mesh%tstate = rk4a(irk)*mesh%tstate + dt*mesh%hstate
      mesh%state = mesh%state + rk4b(irk)*mesh%tstate
    else if (timeIntegrationMethod == 2) then

      tu = rk3a(irk)*tu + dt*hu
      u = u + rk3b(irk)*tu

      mesh%tslip1 = rk3a(irk)*mesh%tslip1 + dt*mesh%sliprate1
      mesh%slip1 = mesh%slip1 + rk3b(irk)*mesh%tslip1
      mesh%tslip2 = rk3a(irk)*mesh%tslip2 + dt*mesh%sliprate2
      mesh%slip2 = mesh%slip2 + rk3b(irk)*mesh%tslip2
      ! rate state
      mesh%tstate = rk3a(irk)*mesh%tstate + dt*mesh%hstate
      mesh%state = mesh%state + rk3b(irk)*mesh%tstate
    else if (timeIntegrationMethod == 3) then

      tu = rk2a(irk)*tu + dt*hu
      u = u + rk2b(irk)*tu

      mesh%tslip1 = rk2a(irk)*mesh%tslip1 + dt*mesh%sliprate1
      mesh%slip1 = mesh%slip1 + rk2b(irk)*mesh%tslip1
      mesh%tslip2 = rk2a(irk)*mesh%tslip2 + dt*mesh%sliprate2
      mesh%slip2 = mesh%slip2 + rk2b(irk)*mesh%tslip2
      ! rate state
      mesh%tstate = rk2a(irk)*mesh%tstate + dt*mesh%hstate
      mesh%state = mesh%state + rk2b(irk)*mesh%tstate

    end if

    mesh%slip = sqrt(mesh%slip1**2+mesh%slip2**2)

    !if (mod(it-1,grdsurf_snap_skip) == 0) then
    !  call grdsurf_io_save(mesh,u,displ,(it-1)/grdsurf_snap_skip+1)
    !end if

    if (use_damp == 1) then
      do i = 1,Nvar
        u(:,i) = u(:,i)*mesh%damp
      end do
    end if

  end do ! end rk

  ! displacement
  displ = displ + dt*(mu+u)*0.5

  if(plasticity==1) &
  call update_plastic(mesh,u)

  call update_info(mesh)

  ! thermal pressurization
  if(thermalpressure==1) &
  call add_TP(mesh)

  if(export_wave ==1 .and. it==export_wave_timestep) then
    call write_wave_vtk(mesh,u)
  end if

  call energy_io_save(mesh,u,it)

end do

call fault_io_end(mesh)
call grdsurf_io_end(mesh)
call recv_io_end(mesh)
call energy_io_end(mesh)

call cpu_time(time_end)

print*, 'cost:',time_end-time_start,' sec'
call MPI_Finalize(ierr)

END PROGRAM
