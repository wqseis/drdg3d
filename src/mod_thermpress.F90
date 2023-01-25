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

module mod_thermpress

  use mod_para, only : RKIND,       &
                       Nfp, Nfaces, &
                       BC_FAULT,    &
                       PI, SQRT2PI
  use mod_types

  implicit none


contains

subroutine add_TP(mesh)
  implicit none

  type(meshvar) :: mesh

  integer :: ie,ief,is,i,it

  real(kind=RKIND) :: Lambda,rho_c,alpha_th,alpha_hy,TP_w
  integer :: TP_n
  real(kind=RKIND) :: temp,pressure
  real(kind=RKIND) :: tauV
  real(kind=RKIND) :: dt
  real(kind=RKIND) :: theta(mesh%TP_n),sigma(mesh%TP_n)

  dt = mesh%deltat

  Lambda = 0.1 ! 0.1 MPa/K 
  rho_c = 2.7  ! 2.7 MJ/m^3/K (MPa/K)
  alpha_th = 1d-6 ! 1e-6 m^2/s
  alpha_hy = 4d-4 ! 4e-4 m^2/s
  TP_w = 0.02     ! 20 mm
  TP_n = mesh%TP_n

  do ief = 1,mesh%nfault_elem
    ie = mesh%fault2wave(ief)
    do is = 1,Nfaces
      if (mesh%bctype(is,ie) >= BC_FAULT) then
        do i = 1,Nfp
          alpha_hy = mesh%TP_hy(i,is,ief)
          do it = 1,1
          tauV = mesh%stress(i,is,ief) * mesh%sliprate(i,is,ief)
          !auV = 10.0*1.0
          theta = mesh%TP_Theta(:,i,is,ief)
          sigma = mesh%TP_Sigma(:,i,is,ief)
          call calc_thermpress(dt/1.0,TP_n,TP_w,alpha_th,alpha_hy,rho_c,Lambda,&
              theta,sigma,tauV,mesh%Dwn,mesh%DFinv,temp,pressure)

          mesh%TP_Theta(:,i,is,ief) = theta
          mesh%TP_Sigma(:,i,is,ief) = sigma

          mesh%TP_T(i,is,ief) = temp
          mesh%TP_P(i,is,ief) = pressure
          end do
        end do
      end if
    end do
  end do

end subroutine

subroutine calc_thermpress(dt,TP_n,TP_w,alpha_th,alpha_hy,rho_c,Lambda,&
    theta,sigma,tauV,Dwn,DFinv,temp,pressure)
  implicit none
  !type(meshvar) :: mesh
  real(kind=RKIND) :: dt
  integer :: TP_n
  real(kind=RKIND) :: TP_w
  real(kind=RKIND) :: alpha_th,alpha_hy,rho_c
  real(kind=RKIND) :: Lambda
  real(kind=RKIND) :: theta(TP_n),sigma(TP_n)
  real(kind=RKIND) :: theta_current(TP_n),sigma_current(TP_n)
  real(kind=RKIND) :: tauV ! shear stress * slip rate
  real(kind=RKIND) :: Dwn(TP_n),DFinv(TP_n)
  real(kind=RKIND) :: temp,pressure
  ! local variables
  real(kind=RKIND) :: Lambda_prime
  real(kind=RKIND) :: tmp(TP_n),x_th(TP_n),x_hy(TP_n)!,x(TP_n)
  real(kind=RKIND) :: omega(TP_n) ! shear heating source
  real(kind=RKIND) :: T,P
  !integer :: i
  real(kind=RKIND) :: temp_0, pressure_0

  intent(in) :: dt,TP_n,TP_w,alpha_th,alpha_hy,rho_c,Lambda,tauV,Dwn,DFinv
  intent(inout) :: theta,sigma
  intent(out) :: temp,pressure

  Lambda_prime = Lambda*alpha_th/(alpha_hy-alpha_th)
  tmp = (Dwn/TP_w)**2
  x_th = alpha_th*dt*tmp
  x_hy = alpha_hy*dt*tmp

  ! 1. calculate diffusion at previous timestep
  ! temperature
  theta_current = theta*exp(-alpha_th*dt*tmp)
  ! pore pressure + lambda' * temp
  sigma_current = sigma*exp(-alpha_hy*dt*tmp)

  ! 2. add current contribution and get new temperature
  call heat_source(TP_w,alpha_th,dt,Dwn,TP_n,omega)
  theta = theta_current + (tauV/rho_c)*omega
  call heat_source(TP_w,alpha_hy,dt,Dwn,TP_n,omega)
  sigma = sigma_current + ((Lambda+Lambda_prime)*tauV)/rho_c*omega

  ! 3. recover temperature and pressure using iFFT
  T = 0.0; P = 0.0
  ! new contribution
  !do i = 1,TP_n
  !  T = T + (DFinv(i)/TP_w)*theta(i)
  !  P = P + (DFinv(i)/TP_w)*sigma(i)
  !end do
  T = dot_product(DFinv,theta)/TP_w
  P = dot_product(DFinv,sigma)/TP_w

  ! update pore pressure change (sigma = pore pressure + lambda' * temp)
  P = P - Lambda_prime*T

  temp_0 = 0.0
  pressure_0 = 0.0

  ! temperature and pore pressure change
  temp = T + temp_0
  pressure = p + pressure_0

end subroutine

  !call heat_source(TP_w,alpha_th,dt,Dwn,TP_n,omega)
subroutine heat_source(w,alpha,dt,Dwn,TP_n,omega)
  implicit none
  !real(kind=RKIND),parameter :: SQRT2PI = 2.506628274631000
  real(kind=RKIND) :: w ! half width of the shear zone
  real(kind=RKIND) :: alpha,dt
  integer :: TP_n ! number of points of the shear zone (half)
  real(kind=RKIND) :: Dwn(TP_n)
  real(kind=RKIND) :: omega(TP_n)
  ! local variables
  real(kind=RKIND) :: tmp(TP_n)

  intent(in) :: w,alpha,dt,Dwn,TP_n
  intent(out) :: omega

  ! Gaussian shear zone in spectral domain
  tmp = (Dwn/w)**2
  !tmp = alpha*dt*tmp
  ! heat source function in the wavenumber domain
  !omega = (1d0-dexp(-tmp))/tmp
  !omega = omega*dexp(-0.5*Dwn**2)/SQRT2PI*dt
  omega = 1.0/(alpha*tmp*SQRT2PI)*exp(-0.5*(Dwn)**2)*(1.0-exp(-alpha*dt*tmp))
end subroutine

!subroutine init_thermpress(DFinv,Dwn,TP_n)
subroutine init_thermpress(mesh)
  implicit none
  !integer,intent(in) :: TP_n
  integer :: TP_n
  !real(kind=RKIND),intent(inout) :: DFinv(TP_n),Dwn(TP_n)
  type(meshvar) :: mesh

  real(kind=RKIND) :: TP_log_dz, TP_max_wavenumber
  real(kind=RKIND) :: TP_grid

  integer :: j

  TP_n = 60

  mesh%TP_n = TP_n

  allocate(mesh%DFinv(mesh%TP_n))
  allocate(mesh%Dwn  (mesh%TP_n))
  allocate(mesh%TP_Theta(mesh%TP_n,Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%TP_Sigma(mesh%TP_n,Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%TP_T(Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%TP_P(Nfp,Nfaces,mesh%nfault_elem))
  !allocate(mesh%TP_hy(Nfp,Nfaces,mesh%nfault_elem))

  mesh%TP_Theta(:,:,:,:) = 0.0
  mesh%TP_Sigma(:,:,:,:) = 0.0
  mesh%TP_T(:,:,:) = 0.0
  mesh%TP_P(:,:,:) = 0.0


  ! Noda and Lapusta (2010)
  TP_log_dz = 0.3d0
  TP_max_wavenumber = 10.0d0

  do j = 1,TP_n
    TP_grid = TP_max_wavenumber*dexp(-TP_log_dz*dble(TP_n-j))
    mesh%Dwn(j) = TP_grid
    if (j == 1) then
      mesh%DFinv(j) = dsqrt(2d0/PI)*TP_grid*(1d0+TP_log_dz*0.5)
      !mesh%DFinv(j) = sqrt(2.0/PI)*TP_grid*(1.0+TP_log_dz)
    elseif (j == TP_n) then
      mesh%DFinv(j) = dsqrt(2d0/PI)*TP_grid*TP_log_dz*0.5
    else
      mesh%DFinv(j) = dsqrt(2d0/PI)*TP_grid*TP_log_dz
    end if
  end do
end subroutine

end module
