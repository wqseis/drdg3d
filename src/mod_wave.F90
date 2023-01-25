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

module mod_wave

  use mod_para,   only : RKIND, PI,                   &
                         Np, Nfp, Nfaces, Nvar,       &
                         initial_condition_wave,      &
                         src_loc,                     &
                         src_gaussian_width,          &
                         src_m0,                      &
                         src_mxx,                     &
                         src_myy,                     &
                         src_mzz,                     &
                         src_myz,                     &
                         src_mxz,                     &
                         src_mxy,                     &
                         problem,                     &
                         mesh_dir,                    &
                         thermalpressure,             &
                         smooth_load,                 &
                         smooth_load_time,            &
                         friction_law,                &
                         RS_f0,                       &
                         RS_V0,                       &
                         RS_fw,                       &
                         flux_method,                 &
                         BC_IN, BC_FREE, BC_FAULT,    &
                         nucleate_y0,                 &
                         nucleate_z0,                 &
                         nucleate_Vrup,               &
                         nucleate_rcrit,              &
                         TimeForcedRup,               &
                         export_wave_component
  use mod_types,  only : meshvar
  use mod_rotate, only : rotate_xyz2nml,              &
                         rotate_nml2xyz,              &
                         rotate_u,                    &
                         rotation_matrix_strain,      &
                         rotation_matrix_strain_inv,  &
                         rotation_matrix_velocity,    &
                         rotation_matrix_velocity_inv
  use mod_eqns,   only : strain2stress
  use mod_solve,  only : Regula_Falsi
  use mod_vtk,    only : writeVtkTetraMeshRealdata
  use mod_eqns,   only : Flux1,Flux2,Flux3,           &
                         extract_traction_velocity,   &
                         riemannSolver_continuous,    &
                         generate_fluctuations

  use mod_math,   only : mxm

  implicit none

contains

subroutine RHS(mesh,u,qi,ru)
  implicit none
  type(meshvar) :: mesh
  real(kind=RKIND),dimension(:,:) :: u,ru ! (Np*Nelem,Nvar)
  real(kind=rkind),intent(in) :: qi(:,:,:,:) ! (...
  real(kind=RKIND),dimension(Np,Nvar) :: Ue!,F1,F2,F3 ! (Np,Nvar)
  !real(kind=RKIND),dimension(Np) :: dF1dr,dF1ds,dF1dt
  !real(kind=RKIND),dimension(Np) :: dF2dr,dF2ds,dF2dt
  !real(kind=RKIND),dimension(Np) :: dF3dr,dF3ds,dF3dt
  real(kind=RKIND),dimension(Nfp*Nfaces,Nvar) :: fluxs
  real(kind=RKIND) :: cp,cs,rho,rrho
  real(kind=RKIND),dimension(Np) :: vx,vy,vz
  real(kind=RKIND) :: exx,eyy,ezz,eyz,exz,exy
  real(kind=RKIND),dimension(Np) :: sxx,syy,szz,syz,sxz,sxy
  !real(kind=RKIND),dimension(Np) :: Trx,Try,Trz
  !real(kind=RKIND),dimension(Np) :: Tsx,Tsy,Tsz
  !real(kind=RKIND),dimension(Np) :: Ttx,Tty,Ttz
  !real(kind=RKIND),dimension(Np) :: DrTrx,DrTry,DrTrz
  !real(kind=RKIND),dimension(Np) :: DsTsx,DsTsy,DsTsz
  !real(kind=RKIND),dimension(Np) :: DtTtx,DtTty,DtTtz
  real(kind=RKIND),dimension(Np) :: DrSxx,DrSyy,DrSzz,DrSyz,DrSxz,DrSxy,DrVx,DrVy,DrVz
  real(kind=RKIND),dimension(Np) :: DsSxx,DsSyy,DsSzz,DsSyz,DsSxz,DsSxy,DsVx,DsVy,DsVz
  real(kind=RKIND),dimension(Np) :: DtSxx,DtSyy,DtSzz,DtSyz,DtSxz,DtSxy,DtVx,DtVy,DtVz
  real(kind=RKIND),dimension(Np,Nvar) :: ru1
  real(kind=RKIND),dimension(Np) :: rx,ry,rz
  real(kind=RKIND),dimension(Np) :: sx,sy,sz
  real(kind=RKIND),dimension(Np) :: tx,ty,tz
  real(kind=RKIND) :: dt,depth
  !integer :: iv(Np)
  integer :: Nelem,ie,i,j,k,k1,k2

  !allocate(Ue(Np,Nvar))
  !allocate(F1(Np,Nvar))
  !allocate(F2(Np,Nvar))
  !allocate(F3(Np,Nvar))
  !allocate(dF1dr(Np))
  !allocate(dF1ds(Np))
  !allocate(dF1dt(Np))
  !allocate(dF2dr(Np))
  !allocate(dF2ds(Np))
  !allocate(dF2dt(Np))
  !allocate(dF3dr(Np))
  !allocate(dF3ds(Np))
  !allocate(dF3dt(Np))
  !allocate(fluxs(Nfp*Nfaces,Nvar))
  !allocate(flux(Nvar))

  Nelem = mesh%Nelem

  dt = mesh%deltat

  do ie = 1,Nelem
    !do j = 1,Np
    !  iv(j) = j+(ie-1)*Np
    !end do
    cp = mesh%vp(ie)
    cs = mesh%vs(ie)
    rho = mesh%rho(ie)
    rrho = 1d0/rho

    do j = 1,Np
      k = j+(ie-1)*Np
      !Ue(j,:) = u(j+(ie-1)*Np,:)
      Ue(j,:) = u(k,:)

      rx(j) = mesh%rx(k)
      ry(j) = mesh%ry(k)
      rz(j) = mesh%rz(k)

      sx(j) = mesh%sx(k)
      sy(j) = mesh%sy(k)
      sz(j) = mesh%sz(k)

      tx(j) = mesh%tx(k)
      ty(j) = mesh%ty(k)
      tz(j) = mesh%tz(k)
    end do

    ! U=(rho*Vx,rho*Vy,rho*Vz,Exx,Eyy,Ezz,Eyz,Exz,Exy)
    ! F1=(Sxx,Sxy,Sxz,Vx,0,0,0,Vz,Vy)
    ! F2=(Sxy,Syy,Syz,0,Vy,0,Vz,0,Vx)
    ! F3=(Sxz,Syz,Szz,0,0,Vz,Vy,Vx,0)
    do i = 1,Np

      depth = -mesh%vz(i+(ie-1)*Np)

      vx(i)  = Ue(i,1)*rrho
      vy(i)  = Ue(i,2)*rrho
      vz(i)  = Ue(i,3)*rrho
      exx = Ue(i,4)
      eyy = Ue(i,5)
      ezz = Ue(i,6)
      eyz = Ue(i,7)
      exz = Ue(i,8)
      exy = Ue(i,9)
      call strain2stress(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs, &
                         sxx(i),syy(i),szz(i),syz(i),sxz(i),sxy(i))

      !if(mesh%irk==mesh%nrk) call Return_Map(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs,depth,dt,&
      !                                       sxx,syy,szz,syz,sxz,sxy)

      !F1(i,:) = 0; F2(i,:) = 0; F3(i,:) = 0
      !F1(i,1) = sxx; F1(i,2) = sxy; F1(i,3) = sxz;
      !F2(i,1) = sxy; F2(i,2) = syy; F2(i,3) = syz;
      !F3(i,1) = sxz; F3(i,2) = syz; F3(i,3) = szz;
      !
      !F1(i,4) = Vx
      !F1(i,8) = Vz
      !F1(i,9) = Vy

      !F2(i,5) = Vy
      !F2(i,7) = Vz
      !F2(i,9) = Vx

      !F3(i,6) = Vz
      !F3(i,7) = Vy
      !F3(i,8) = Vx

      !call Flux1(Ue(i,:),rho,cp,cs,F1(i,:))
      !call Flux2(Ue(i,:),rho,cp,cs,F2(i,:))
      !call Flux3(Ue(i,:),rho,cp,cs,F3(i,:))
    end do

    !do i = 1,Nvar
    !  dF1dr = matmul(mesh%Dr,F1(:,i))
    !  dF1ds = matmul(mesh%Ds,F1(:,i))
    !  dF1dt = matmul(mesh%Dt,F1(:,i))
    !  dF2dr = matmul(mesh%Dr,F2(:,i))
    !  dF2ds = matmul(mesh%Ds,F2(:,i))
    !  dF2dt = matmul(mesh%Dt,F2(:,i))
    !  dF3dr = matmul(mesh%Dr,F3(:,i))
    !  dF3ds = matmul(mesh%Ds,F3(:,i))
    !  dF3dt = matmul(mesh%Dt,F3(:,i))

    !  k1 = 1+(ie-1)*Np
    !  k2 = Np+(ie-1)*Np
    !  ru(k1:k2,i) = &
    !    mesh%rx(k1:k2)*dF1dr+mesh%sx(k1:k2)*dF1ds+mesh%tx(k1:k2)*dF1dt + &
    !    mesh%ry(k1:k2)*dF2dr+mesh%sy(k1:k2)*dF2ds+mesh%ty(k1:k2)*dF2dt + &
    !    mesh%rz(k1:k2)*dF3dr+mesh%sz(k1:k2)*dF3ds+mesh%tz(k1:k2)*dF3dt
    !end do

    !@ TrX = Sxx*rx+Sxy*ry+Sxz*rz
    !@ TrY = Sxy*rx+Sxy*ry+Sxz*rz
    !@ TrZ = Sxz*rx+Sxy*ry+Sxz*rz

    !@ TsX = Sxx*sx+Sxy*sy+Sxz*sz
    !@ TsY = Sxy*sx+Sxy*sy+Sxz*sz
    !@ TsZ = Sxz*sx+Sxy*sy+Sxz*sz

    !@ TtX = Sxx*tx+Sxy*ty+Sxz*tz
    !@ TtY = Sxy*tx+Sxy*ty+Sxz*tz
    !@ TtZ = Sxz*tx+Sxy*ty+Sxz*tz

    !@ DrTrx = matmul(mesh%Dr,Trx)
    !@ DrTry = matmul(mesh%Dr,Try)
    !@ DrTrz = matmul(mesh%Dr,Trz)
    !@ DsTsx = matmul(mesh%Ds,Tsx)
    !@ DsTsy = matmul(mesh%Ds,Tsy)
    !@ DsTsz = matmul(mesh%Ds,Tsz)
    !@ DtTtx = matmul(mesh%Dt,Ttx)
    !@ DtTty = matmul(mesh%Dt,Tty)
    !@ DtTtz = matmul(mesh%Dt,Ttz)

    DrSxx = matmul(mesh%Dr,Sxx)
    DrSyy = matmul(mesh%Dr,Syy)
    DrSzz = matmul(mesh%Dr,Szz)
    DrSyz = matmul(mesh%Dr,Syz)
    DrSxz = matmul(mesh%Dr,Sxz)
    DrSxy = matmul(mesh%Dr,Sxy)
    DrVx  = matmul(mesh%Dr,Vx )
    DrVy  = matmul(mesh%Dr,Vy )
    DrVz  = matmul(mesh%Dr,Vz )

    DsSxx = matmul(mesh%Ds,Sxx)
    DsSyy = matmul(mesh%Ds,Syy)
    DsSzz = matmul(mesh%Ds,Szz)
    DsSyz = matmul(mesh%Ds,Syz)
    DsSxz = matmul(mesh%Ds,Sxz)
    DsSxy = matmul(mesh%Ds,Sxy)
    DsVx  = matmul(mesh%Ds,Vx )
    DsVy  = matmul(mesh%Ds,Vy )
    DsVz  = matmul(mesh%Ds,Vz )

    DtSxx = matmul(mesh%Dt,Sxx)
    DtSyy = matmul(mesh%Dt,Syy)
    DtSzz = matmul(mesh%Dt,Szz)
    DtSyz = matmul(mesh%Dt,Syz)
    DtSxz = matmul(mesh%Dt,Sxz)
    DtSxy = matmul(mesh%Dt,Sxy)
    DtVx  = matmul(mesh%Dt,Vx )
    DtVy  = matmul(mesh%Dt,Vy )
    DtVz  = matmul(mesh%Dt,Vz )

    ru1(:,1) = &
      DrSxx*rx+DsSxx*sx+DtSxx*tx+&
      DrSxy*ry+DsSxy*sy+DtSxy*ty+&
      DrSxz*rz+DsSxz*sz+DtSxz*tz
    ru1(:,2) = &
      DrSxy*rx+DsSxy*sx+DtSxy*tx+&
      DrSyy*ry+DsSyy*sy+DtSyy*ty+&
      DrSyz*rz+DsSyz*sz+DtSyz*tz
    ru1(:,3) = &
      DrSxz*rx+DsSxz*sx+DtSxz*tx+&
      DrSyz*ry+DsSyz*sy+DtSyz*ty+&
      DrSzz*rz+DsSzz*sz+DtSzz*tz
    !ru1(:,1) = DrTrx+DsTsx+DtTtx
    !ru1(:,2) = DrTry+DsTsy+DtTty
    !ru1(:,3) = DrTrz+DsTsz+DtTtz
    ru1(:,4) = DrVx*rx+DsVx*sx+DtVx*tx
    ru1(:,5) = DrVy*ry+DsVy*sy+DtVy*ty
    ru1(:,6) = DrVz*rz+DsVz*sz+DtVz*tz
    ru1(:,7) = &
         DrVz*ry+DsVz*sy+DtVz*ty + &
         DrVy*rz+DsVy*sz+DtVy*tz
    ru1(:,8) = &
         DrVz*rx+DsVz*sx+DtVz*tx + &
         DrVx*rz+DsVx*sz+DtVx*tz
    ru1(:,9) = &
         DrVy*rx+DsVy*sx+DtVy*tx + &
         DrVx*ry+DsVx*sy+DtVx*ty

    do i = 1,Nvar
      do j = 1,Np
        k = j+(ie-1)*Np
        ru(k,i) = ru1(j,i)
      end do
    end do

  !end do ! element


    ! compute numerical flux
    !do is = 1,Nfaces
    !  do i = 1,Nfp
    !    call get_flux(mesh,u,i,is,ie,qi,rho,cp,cs,flux)
    !    fluxs(i+(is-1)*Nfp,:) = flux
    !  end do
    !end do
    call get_flux(mesh,u,ie,qi,fluxs)
    !fluxs=0.
  !call get_flux(mesh,u,qi,rho,cp,cs)

  !do ie = 1,mesh%Nelem
    !j=(ie-1)*Nfp*Nfaces+1
    !k=ie*Nfp*Nfaces
    !fluxs = fluxes(j:k,:)
    !do i = 1,Nvar
    !  k1 = 1+(ie-1)*Np
    !  k2 = Np+(ie-1)*Np
    !  !ru(k1:k2,i) = ru(k1:k2,i) + matmul(mesh%LIFT,mesh%Fscale(:,ie)*fluxs(:,i))
    !  ru(k1:k2,i) = ru(k1:k2,i) + matmul(mesh%LIFT,fluxs(:,i))
    !end do
    !do i = 1,Nvar
      k1 = (ie-1)*Np+1
      k2 = ie*Np
      !ru(k1:k2,i) = ru(k1:k2,i) + matmul(mesh%LIFT,mesh%Fscale(:,ie)*fluxs(:,i))
      !call mxm(mesh%LIFT,fluxs,fluxs1,Np,Nfp*Nfaces,Nvar)
      ru(k1:k2,:) = ru(k1:k2,:) + matmul(mesh%LIFT,fluxs)
      !ru(k1:k2,:) = ru(k1:k2,:) + fluxs1

    !end do

  end do ! element

  !deallocate(Ue)
  !deallocate(F1,F2,F3)
  !deallocate(dF1dr,dF1ds,dF1dt)
  !deallocate(dF2dr,dF2ds,dF2dt)
  !deallocate(dF3dr,dF3ds,dF3dt)
  !deallocate(fluxs,flux)
end subroutine

!subroutine get_flux(mesh,u,i,is,ie,qi,rho,cp,cs,fstar)
subroutine get_flux(mesh,u,ie,qi,fluxs)
!subroutine get_flux(mesh,u,qi,rho,cp,cs)
  implicit none
  type(meshvar) :: mesh
  real(kind=RKIND),dimension(:,:) :: u
  real(kind=RKIND),dimension(:,:,:,:) :: qi
  integer :: i,j,is,ie,ief
  integer :: neigh,direction
  real(kind=RKIND) :: rho,rrho,cp,cs
  real(kind=RKIND) :: rho_out,rrho_out,cp_out,cs_out
  real(kind=RKIND),dimension(Nfp*Nfaces,Nvar) :: fluxs
  real(kind=RKIND),dimension(Nvar) :: uL,uR,fstar!,FL,FR!,du,dF,dFx,dFy,dFz
  real(kind=RKIND),dimension(3) :: n_n,n_m,n_l
  !real(kind=RKIND),dimension(3,3) :: matTv,invTv
  !real(kind=RKIND),dimension(6,6) :: matTs,invTs
  real(kind=RKIND) :: Zp,Zs,Zs_in,Zs_out,Zp_out,eta
  !real(kind=RKIND) :: alpha(3)
  real(kind=RKIND) :: tau0n,tau0m,tau0l,vv1,vv2,vel
  real(kind=RKIND) :: tau_lock,tau_lock_1,tau_lock_2,tau_str,tau_n
  !real(kind=RKIND) :: V_ref,L_ref,coef,tau_n_eff,strength_exp,sigma,dt
  real(kind=RKIND) :: sigma,dt
  real(kind=RKIND) :: fault_mu,mu_s,mu_d,Dc,C0,sliprate,slip
  integer :: flipped_index(Nfp)
  integer :: mpi_e,mpi_n
  !real(kind=RKIND) :: exx,eyy,ezz,eyz,exz,exy
  !real(kind=RKIND) :: sxx,syy,szz,syz,sxz,sxy
  real(kind=RKIND) :: Tx,Ty,Tz!,Tn,Tm,Tl
  real(kind=RKIND) :: Vx,Vy,Vz!,Vn,Vm,Vl
  real(kind=RKIND) :: depth
  ! rate state
  real(kind=RKIND) :: a,b,Vw,psi,Phi,sigma_n,V,V0,f0,L0,fv,flv,fss,fw,psiss
  real(kind=RKIND) :: Gt
  real(kind=RKIND) :: vn_p,vm_p,vl_p,Tn_p,Tm_p,Tl_p
  real(kind=RKIND) :: vn_m,vm_m,vl_m,Tn_m,Tm_m,Tl_m
  real(kind=RKIND) :: vn_hat_p,vm_hat_p,vl_hat_p,Tn_hat_p,Tm_hat_p,Tl_hat_p
  real(kind=RKIND) :: vn_hat_m,vm_hat_m,vl_hat_m,Tn_hat_m,Tm_hat_m,Tl_hat_m
  real(kind=RKIND) :: FLn,FLm,FLl
  !real(kind=RKIND) :: FRn,FRm,FRl
  real(kind=RKIND) :: FL_n,FL_m,FL_l
  !real(kind=RKIND) :: FR_n,FR_m,FR_l
  real(kind=RKIND) :: FLx,FLy,FLz
  !real(kind=RKIND) :: FRx,FRy,FRz
  real(kind=RKIND) :: FL_x,FL_y,FL_z
  !real(kind=RKIND) :: FR_x,FR_y,FR_z
  real(kind=RKIND) :: dVx,dVy,dVz
  real(kind=RKIND) :: dTx,dTy,dTz



!!#if defined (TPV24) || defined (TPV26) || defined (TPV29) || defined (WENCHUAN) || defined (TW)
!!!#if defined (TPV24) || defined (TPV26) || defined (TPV29) || defined (TW)
  !real(kind=RKIND) :: SW_T, SW_t0, SW_f1, SW_f2, dist, SW_rcrit
  real(kind=RKIND) :: x,y,z,cur_time!,Pf
  !SW_rcrit = sqrt(10*10/PI)
  !SW_rcrit = 4.0
  !SW_t0 = 0.5
  !SW_t0 = 0.005
!!!#endif

  !do ie = 1,mesh%Nelem
  ief = mesh%wave2fault(ie) ! =0 if not fault

  rho=mesh%rho(ie)
  cp=mesh%vp(ie)
  cs=mesh%vs(ie)
  rrho = 1d0/rho

  Zp = rho*cp
  Zs = rho*cs

  ! will overwrite
  rho_out = rho
  rrho_out = rrho
  cp_out = cp
  cs_out = cs

  dt = mesh%deltat

  do is = 1,Nfaces
  do i = 1,Nfp
  depth = -mesh%vz(mesh%vmapM(i,is,ie))

  uL = u(mesh%vmapM(i,is,ie),:)

  neigh = mesh%neigh(is,ie)
  direction = mesh%direction(is,ie)
  !if(neigh>0) uR=u(mesh%vmapP(i,is,ie),:)
  if(neigh==ie) uR=0 ! absorbing
  !if(neigh>0 .and. neigh .ne. ie) then
  if(neigh>0) then
    uR=u(mesh%vmapP(i,is,ie),:)
    rho_out = mesh%rho(neigh)
    cp_out = mesh%vp(neigh)
    cs_out = mesh%vs(neigh)
  elseif (neigh==0) then
    uR(:) = 0
  elseif (neigh==-1) then
    !print*,'mpi'
    mpi_e = mesh%mpi_ibool(is,ie)
    mpi_n = mesh%mpi_interface(4,is,ie)
    !print*,'mpi_e=',mpi_e
    !print*,'mpi_n=',mpi_n
    flipped_index = mesh%flipped_index(:,direction)
    uR(:) = qi(flipped_index(i),1:9,mpi_e,mpi_n)  ! mpi
    !Loop over all faces that sit on an mpi-interface for that rank/processor
    !do i = 1, size(mesh%mpi_vp(:,1))
    do j = 1, mesh%pinterfaces
      !Check if the rank and the element in that rank are in the mpi-impedance array (they should)
      if ((abs(mesh%mpi_vp(j,2) - mesh%mpi_interface(1,is,ie)) < epsilon(mesh%mpi_vp(j,2))) .and. &
          (abs(mesh%mpi_vp(j,3) - mesh%mpi_interface(2,is,ie)) < epsilon(mesh%mpi_vp(j,3)))) then
      !if ( int(mesh%mpi_vp(i,2)) .eq. mesh%mpi_interface(1,is,ie) .and. &
      !     int(mesh%mpi_vp(i,3)) .eq. mesh%mpi_interface(2,is,ie) ) then
        !associate corresponding impedance value to zout
        rho_out = mesh%mpi_rho(j,1)
        cp_out = mesh%mpi_vp(j,1)
        cs_out = mesh%mpi_vs(j,1)
        rrho_out = 1d0/rho_out
      end if
    end do
  else
    uR(:) = 1e38
  end if
  if(neigh==ie) uR=0 ! absorbing

  if (mesh%bctype(is,ie) == BC_FREE) then ! free surface
      uR(1:3) = uL(1:3)
      uR(4:9) = -uL(4:9)
      !uR(:) = 0
  end if

  !if (mesh%bctype(is,ie) >= BC_FAULT) then
  !    uR = uL
  !    uR(1:3) = -uL(1:3)
  !end if

  Zp_out = rho_out * cp_out
  Zs_out = rho_out * cs_out

  j = i+(is-1)*Nfp
  n_n = (/mesh%nx(j,ie),mesh%ny(j,ie),mesh%nz(j,ie)/)
  n_m = (/mesh%mx(j,ie),mesh%my(j,ie),mesh%mz(j,ie)/)
  n_l = (/mesh%lx(j,ie),mesh%ly(j,ie),mesh%lz(j,ie)/)

  !n_m = (/1,0,0/)
  !n_n = (/0,1,0/)
  !n_l = (/0,0,1/)
  !print*,n_n
  !if(abs(norm3(n_n)-1.0)>1e-3) then
  !print*,n_n
  !end if

!  du = uR-uL
!  call Flux1(du,rho,cp,cs,dFx)
!  call Flux2(du,rho,cp,cs,dFy)
!  call Flux3(du,rho,cp,cs,dFz)
!  fstar = dFx*n_n(1)+dFy*n_n(2)+dFz*n_n(3)
!  fstar = fstar * 0.5
  !return

  !call rotation_matrix_velocity    (n_n,n_m,n_l,matTv)
  !call rotation_matrix_strain      (n_n,n_m,n_l,matTs)
  !call rotation_matrix_velocity_inv(n_n,n_m,n_l,invTv)
  !call rotation_matrix_strain_inv  (n_n,n_m,n_l,invTs)

  !if (.false.) then
  !call rotate_u(matTv,matTs,uL)
  !call rotate_u(matTv,matTs,uR)

  !call Flux1(uL,rho,cp,cs,FL)
  !call Flux1(uR,rho,cp,cs,FR)
  !end if

  !if (.false.) then
  ! ============= - side ==============================
  !Vx  = uL(1)*rrho
  !Vy  = uL(2)*rrho
  !Vz  = uL(3)*rrho
  !exx = uL(4)
  !eyy = uL(5)
  !ezz = uL(6)
  !eyz = uL(7)
  !exz = uL(8)
  !exy = uL(9)
  !call strain2stress(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs, &
  !                   sxx,syy,szz,syz,sxz,sxy)
  !!if(mesh%irk==mesh%nrk) call Return_Map(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs,depth,dt,&
  !!                                       sxx,syy,szz,syz,sxz,sxy)
  !Tx = sxx*n_n(1)+sxy*n_n(2)+sxz*n_n(3)
  !Ty = sxy*n_n(1)+syy*n_n(2)+syz*n_n(3)
  !Tz = sxz*n_n(1)+syz*n_n(2)+szz*n_n(3)

  call extract_traction_velocity(uL,n_n,rho,cp,cs,Vx,Vy,Vz,Tx,Ty,Tz)

  dVx = -Vx
  dVy = -Vy
  dVz = -Vz
  dTx = -Tx
  dTy = -Ty
  dTz = -Tz

  ! rotate to local
  call rotate_xyz2nml(n_n,n_m,n_l,Tx,Ty,Tz,Tn_m,Tm_m,Tl_m)
  call rotate_xyz2nml(n_n,n_m,n_l,Vx,Vy,Vz,Vn_m,Vm_m,Vl_m)
  !FL = 0
  !FL(1) = Tn_m
  !FL(2) = Tm_m
  !FL(3) = Tl_m
  !FL(4) = Vn_m
  !FL(8) = Vl_m
  !FL(9) = Vm_m

  !! ============= + side ==============================
  !Vx  = uR(1)*rrho_out
  !Vy  = uR(2)*rrho_out
  !Vz  = uR(3)*rrho_out
  !exx = uR(4)
  !eyy = uR(5)
  !ezz = uR(6)
  !eyz = uR(7)
  !exz = uR(8)
  !exy = uR(9)
  !call strain2stress(exx,eyy,ezz,eyz,exz,exy,rho_out,cp_out,cs_out, &
  !                   sxx,syy,szz,syz,sxz,sxy)
  !!if(mesh%irk==mesh%nrk) call Return_Map(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs,depth,dt,&
  !!                                       sxx,syy,szz,syz,sxz,sxy)
  !Tx = sxx*n_n(1)+sxy*n_n(2)+sxz*n_n(3)
  !Ty = sxy*n_n(1)+syy*n_n(2)+syz*n_n(3)
  !Tz = sxz*n_n(1)+syz*n_n(2)+szz*n_n(3)

  call extract_traction_velocity(uR,n_n,rho_out,cp_out,cs_out,Vx,Vy,Vz,Tx,Ty,Tz)

  dVx = dVx + Vx
  dVy = dVy + Vy
  dVz = dVz + Vz
  dTx = dTx + Tx
  dTy = dTy + Ty
  dTz = dTz + Tz

  if (flux_method == 1 .and. &
      mesh%bctype(is,ie) == BC_IN .and. &
      mesh%fluxtype(is,ie) == 1) then
    !fstar = (fR -fL)*0.5
    fstar(1) = dTx
    fstar(2) = dTy
    fstar(3) = dTz
    fstar(4) = n_n(1)*dVx
    fstar(5) = n_n(2)*dVy
    fstar(6) = n_n(3)*dVz
    fstar(7) = n_n(3)*dVy+n_n(2)*dVz
    fstar(8) = n_n(3)*dVx+n_n(1)*dVz
    fstar(9) = n_n(2)*dVx+n_n(1)*dVy
    fstar = 0.5*fstar
    goto 100
  end if


  ! rotate to local
  call rotate_xyz2nml(n_n,n_m,n_l,Tx,Ty,Tz,Tn_p,Tm_p,Tl_p)
  call rotate_xyz2nml(n_n,n_m,n_l,Vx,Vy,Vz,Vn_p,Vm_p,Vl_p)
  !FR = 0
  !FR(1) = Tn_p
  !FR(2) = Tm_p
  !FR(3) = Tl_p
  !FR(4) = Vn_p
  !FR(8) = Vl_p
  !FR(9) = Vm_p

  call riemannSolver_continuous(vn_p,vn_m,Tn_p,Tn_m,Zp_out,Zp,vn_hat_p,vn_hat_m,Tn_hat_p,Tn_hat_m)
  call riemannSolver_continuous(vm_p,vm_m,Tm_p,Tm_m,Zs_out,Zs,vm_hat_p,vm_hat_m,Tm_hat_p,Tm_hat_m)
  call riemannSolver_continuous(vl_p,vl_m,Tl_p,Tl_m,Zs_out,Zs,vl_hat_p,vl_hat_m,Tl_hat_p,Tl_hat_m)


  !end if
  !FL = 0
  !FL(1) = sxx
  !FL(2) = sxy
  !FL(3) = sxz
  !FL(4) = uL(1)/rho
  !FL(8) = uL(3)/rho
  !FL(9) = uL(2)/rho

  !exx = uR(4)
  !eyy = uR(5)
  !ezz = uR(6)
  !eyz = uR(7)
  !exz = uR(8)
  !exy = uR(9)
  !call strain2stress(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs, &
  !                   sxx,syy,szz,syz,sxz,sxy)
  !FR = 0
  !FR(1) = sxx
  !FR(2) = sxy
  !FR(3) = sxz
  !FR(4) = uR(1)/rho
  !FR(8) = uR(3)/rho
  !FR(9) = uR(2)/rho


  !fstar = 0.5*(fR-fL) ! centering flux
  ! F=(Sxx,Sxy,Sxz,Vx,0,0,0,Vz,Vy)
  !fstar(1) = 0.5*(fR(1)-fL(1)) + 0.5*(fR(4)-fL(4))*Zp
  !fstar(2) = 0.5*(fR(2)-fL(2)) + 0.5*(fR(9)-fL(9))*Zs
  !fstar(3) = 0.5*(fR(3)-fL(3)) + 0.5*(fR(8)-fL(8))*Zs
  !fstar(4) = 0.5*(fR(4)-fL(4)) + 0.5*(fR(1)-fL(1))/Zp
  !fstar(5) = 0
  !fstar(6) = 0
  !fstar(7) = 0
  !fstar(8) = 0.5*(fR(8)-fL(8)) + 0.5*(fR(3)-fL(3))/Zs
  !fstar(9) = 0.5*(fR(9)-fL(9)) + 0.5*(fR(2)-fL(2))/Zs

  !dF = FR-FL
  !alpha(1) = (dF(1)+Zp_out*dF(4))/(Zp+Zp_out)
  !alpha(2) = (dF(2)+Zs_out*dF(9))/(Zs+Zs_out)
  !alpha(3) = (dF(3)+Zs_out*dF(8))/(Zs+Zs_out)

  !fstar(1) = alpha(1)*Zp
  !fstar(2) = alpha(2)*Zs
  !fstar(3) = alpha(3)*Zs
  !fstar(4) = alpha(1)
  !fstar(5) = 0
  !fstar(6) = 0
  !fstar(7) = 0
  !fstar(8) = alpha(3)
  !fstar(9) = alpha(2)

!!!#ifdef WITH_RUPTURE
!  if (flux_method == 1 .and. &
!      mesh%bctype(is,ie) == BC_IN .and. &
!      mesh%fluxtype(is,ie)==1) then
!    !fstar = (fR -fL)*0.5
!  end if

  if (mesh%bctype(is,ie) >= BC_FAULT) then

    zs_in  = rho*cs
    zs_out = rho_out*cs_out
    eta = zs_in*zs_out/(zs_in+zs_out)

    if (smooth_load == 1) then
      Gt = Gt_func(mesh%current_time,smooth_load_time)
    else
      Gt = 1.0
    end if

    Tau0n    = mesh%Tau0n(i,is,ief) + Gt * mesh%dtau0n(i,is,ief)
    Tau0m    = mesh%Tau0m(i,is,ief) + Gt * mesh%dtau0m(i,is,ief)
    Tau0l    = mesh%Tau0l(i,is,ief) + Gt * mesh%dtau0l(i,is,ief)

    slip     = mesh%slip    (i,is,ief)
    sliprate = mesh%sliprate(i,is,ief)
    sigma    = mesh%sigma   (i,is,ief)
    mu_s     = mesh%mu_s    (i,is,ief)
    mu_d     = mesh%mu_d    (i,is,ief)
    Dc       = mesh%Dc      (i,is,ief)
    C0       = mesh%C0      (i,is,ief)

    ! f (sxx,sxy,vx,0,vy)
    ! F=(Sxx,Sxy,Sxz,Vx,0,0,0,Vz,Vy)
    !Tau_lock_1 = fL(2) + fstar(2) + Tau0m
    !Tau_lock_2 = fL(3) + fstar(3) + Tau0l
    Tau_lock_1 = Tm_hat_m + Tau0m
    Tau_lock_2 = Tl_hat_m + Tau0l
    Tau_lock = sqrt(Tau_lock_1**2 + Tau_lock_2**2)
    !Tau_n = fL(1) + fstar(1) + Tau0n ! sxx
    Tau_n = Tn_hat_m + Tau0n
    !Tn_m = sign(1d0,Tau_lock) * Tn_m
    !if (Tau_n > 0.0) Tau_n = 0.0
!!!#ifdef TPV24
!!!      Tau_n = Tau_n + Pf
!!!#endif

    if (thermalpressure == 1) then
      ! thermal pressurization
      Tau_n = Tau_n + mesh%TP_P(i,is,ief)
    end if

#ifdef TPV6
!!!if (.false.) then
!!!      ! for tpv6
!!!      ! Prakash-Clifton
!!!      V_ref = 1.0
!!!      L_ref = 0.2*Dc
!!!      !L_ref = 5*Dc
!!!      dt = mesh%deltat
!!!      coef = dt/L_ref
!!!      Tau_n_eff = Tau_n + exp(-(sliprate+V_ref)*coef)*(mesh%sigma(i,is,ie)-Tau_n)
!!!      !
!!!      if(sliprate > 0.01) Tau_n = Tau_n_eff
!!!end if
#endif


!!!#if defined (TPV24) || defined (TPV26) || defined (TPV29) || defined (WENCHUAN) || defined(TW)
!!!#if defined (TPV24) || defined (TPV26) || defined (TPV29)

    cur_time = mesh%current_time
    x = mesh%vx(mesh%vmapM(i,is,ie))
    y = mesh%vy(mesh%vmapM(i,is,ie))
    z = mesh%vz(mesh%vmapM(i,is,ie))


    ! slip weakening
    fault_mu = mu_s - (mu_s-mu_d) * min(slip/Dc,1.0)

    if (friction_law == 4) then
      ! time weakening
      call time_weakening(y,z,nucleate_y0,nucleate_z0,nucleate_rcrit,&
          nucleate_Vrup,TimeForcedRup,cur_time,&
          slip,Dc,mu_s,mu_d,fault_mu)
    end if

    !!Pf = 9.8*(-z) ! In MPa
    !!Pf = 0.0
    !!Tau_n = Tau_n + Pf ! effective normal stress
    !dist = sqrt( (y+8.0)**2 + (z+10.0)**2) ! TPV24
    !!dist = sqrt( (y+5.0)**2 + (z+10.0)**2) ! TPV26 or TPV29
    !!dist = sqrt( (y-0.0)**2 + (x-20.0)**2) ! TW dip 12
    !!dist = sqrt( (y+.0)**2 + (z+10.0)**2) ! WENCHUAN
    !if(dist<SW_rcrit) then
    !  SW_T=dist/(0.7*cs)+0.081*SW_rcrit/(0.7*cs)*(1./(1.-(dist/SW_rcrit)**2) -1.0)
    !else
    !  SW_T = 1.0e9
    !endif

    !if(slip<Dc) then
    !  SW_f1 = slip/Dc
    !else
    !  SW_f1 = 1.0
    !endif
    !if(cur_time<SW_T) then
    !  SW_f2 = 0.0
    !elseif(cur_time<SW_T+SW_t0) then
    !  SW_f2 = (cur_time-SW_T)/SW_t0
    !else
    !  SW_f2 = 1.0
    !endif

    !fault_mu = mu_s-(mu_s-mu_d)*max(SW_f1, SW_f2)
!#else
      ! slip weakening
!    fault_mu = mu_s - (mu_s-mu_d) * min(slip/Dc,1.0)
!#endif

    Tau_str = fault_mu * max(0.0,-Tau_n) + C0

    !call prakash_cliff_fric(Tau_str,sigma,sliprate,V_ref,L_ref,fault_mu,dt)

    if (abs(Tau_lock) .gt. abs(Tau_str)) then ! fault is slipping
      ! solve for slip rate
      Vel = (abs(Tau_lock)-abs(Tau_str))/eta
      ! same sign with Tau_lock, parallel condition
      vv1 = Tau_lock_1*Vel/(eta*Vel+Tau_str)
      vv2 = Tau_lock_2*Vel/(eta*Vel+Tau_str)
    else
      vv1 = 0.0
      vv2 = 0.0
    endif

    if (friction_law == 1 .or. &
        friction_law == 2 .or. &
        friction_law == 3) then
      ! rate state
      !Phi = sqrt(phi_1**2 + phi_2**2) ! stress-transfer functional
      Phi = Tau_lock
      sigma_n = max(0.0,-Tau_n)
      !V0 = 1e-6
      !f0 = 0.6
      !fw = 0.2
      V0 = RS_V0
      f0 = RS_f0
      fw = RS_fw
      a = mesh%a(i,is,ief)
      b = mesh%b(i,is,ief)
      !if (friction_law == 3) then
      !  Vw = mesh%Vw(i,is,ief)
      !end if
      L0 = mesh%Dc(i,is,ief)
      psi = mesh%state(i,is,ief)
      !if(ief==1 .and. i==1) print*,L0,psi
      !if(abs(y-0)<1 .and. abs(z-0)<1) then
      !  print*,'phi=',Phi,'sigma=',sigma_n,'psi=',psi,'a=',a,'eta=',eta,'V0=',V0,'V=',V
      !end if
      ! initialize V
      !V = sqrt(v1**2 + v2**2)
      V = mesh%sliprate(i,is,ief)
      if (V > Phi) V = 0.5d0*Phi/eta
      ! solve a nonlinear problem for slip-rate: V
      call Regula_Falsi(V,Phi,eta,sigma_n,psi,V0,a)
      ! compute slip velocities
      fv = sigma_n*a*asinh(0.5d0*V/V0*exp(psi/a))
      vv1 = tau_lock_1*V/(eta*V+fv)
      vv2 = tau_lock_2*V/(eta*V+fv)

      if (friction_law == 1) then
        ! ageing law
        mesh%hstate(i,is,ief) = b*V0/L0*exp(-(psi-f0)/b) - V*b/L0
      end if

      if (friction_law == 2) then
        ! slip law
        mesh%hstate(i,is,ief) = -b*V/L0*(log(V/V0)+(psi-f0)/b)
      end if

      if (friction_law == 3) then
        ! slip law with flash heating
        Vw = mesh%Vw(i,is,ief)
        flv = f0 - (b-a)*log(V/V0)
        fss = fw + (flv-fw)/( (1.0+(V/Vw)**8)**0.125 )
        !fss = flv
        psiss = a*(log(sinh(fss/a)) + log(2.0*(V0/V)) )
        mesh%hstate(i,is,ief) = -V/L0*(psi-psiss)
      end if
    end if

      !print*,vv
      !mesh%sliprate(i,is,ie) = abs(vv)!abs(V_p-V_m)
    mesh%sliprate(i,is,ief) = sqrt(vv1**2+vv2**2) !abs(V_p-V_m)
    if (mesh%sliprate(i,is,ief) > 1e-3 .and. mesh%ruptime(i,is,ief) < 0) then
        mesh%ruptime(i,is,ief) = mesh%current_time
    end if

    ! f (sxx,sxy,vx,0,vy)
    ! F=(Sxx,Sxy,Sxz,Vx,0,0,0,Vz,Vy)
    !fstar(2) = fstar(2) - eta*vv1       ! sxy
    !fstar(3) = fstar(3) - eta*vv2       ! sxz
    !fstar(8) = fstar(8) - eta*vv2/zs_in ! vz
    !fstar(9) = fstar(9) - eta*vv1/zs_in ! vy

    Tm_hat_m = Tm_hat_m - eta*vv1
    Tl_hat_m = Tl_hat_m - eta*vv2       ! sxz
    vl_hat_m = vl_hat_m - eta*vv2/zs_in ! vz
    vm_hat_m = vm_hat_m - eta*vv1/zs_in ! vy

    !fstar(1) = Tau_n - Tau0n - fL(1) 

    !fstar(2) =  - eta*vv1       ! sxy
    !fstar(3) =  - eta*vv2       ! sxz
    !fstar(8) =  - eta*vv2/zs_in ! vz
    !fstar(9) =  - eta*vv1/zs_in ! vy

    !mesh%stress(i,is,ief) = sqrt( &
    !        (fstar(2)+fL(2)+mesh%tau0m(i,is,ief))**2+ &
    !        (fstar(3)+fL(3)+mesh%tau0l(i,is,ief))**2)
    !mesh%sigma(i,is,ief) = mesh%tau0n(i,is,ief) + fstar(1) + fL(1)
    !mesh%stress(i,is,ief) = sqrt( &
    !        (fstar(2)+fL(2)+tau0m)**2+ &
    !        (fstar(3)+fL(3)+tau0l)**2)
    !mesh%sigma(i,is,ief) = tau0n + fstar(1) + fL(1)
    mesh%stress(i,is,ief) = sqrt( &
            (Tm_hat_m+tau0m)**2+ &
            (Tl_hat_m+tau0l)**2)

    !mesh%stress(i,is,ie) = Tau_str
    mesh%sigma(i,is,ief) = Tau_n

    !mesh%sliprate1(i,is,ief) = fL(9)+fstar(9)
    !mesh%sliprate2(i,is,ief) = fL(8)+fstar(8)

    mesh%sliprate1(i,is,ief) = vv1
    mesh%sliprate2(i,is,ief) = vv2

    !! save fault recvs
    !do n = 1,mesh%nrecv
    !  if (ie == mesh%recv_elem(n) .and. is == mesh%recv_face(n)) then
    !    mesh%recv_buffer(i,n,1) = mesh%sliprate(i,is,ief)
    !    mesh%recv_buffer(i,n,2) = mesh%stress(i,is,ief)
    !    mesh%recv_buffer(i,n,3) = mesh%sigma(i,is,ief)
    !    mesh%recv_buffer(i,n,4) = mesh%slip(i,is,ief)
    !  end if
    !end do

  end if ! end of BC_FAULT
!!!#endif

  !call rotate_u(invTv,invTs,fstar)

  !call generate_fluctuations(Zp,Tn_m,Tn_hat_m,vn_m,vn_hat_m,FLn)
  !call generate_fluctuations(Zs,Tm_m,Tm_hat_m,vm_m,vm_hat_m,FLm)
  !call generate_fluctuations(Zs,Tl_m,Tl_hat_m,vl_m,vl_hat_m,FLl)

  FLn = 0.5*((Tn_hat_m-Tn_m)+Zp*(vn_hat_m-vn_m))
  FLm = 0.5*((Tm_hat_m-Tm_m)+Zs*(vm_hat_m-vm_m))
  FLl = 0.5*((Tl_hat_m-Tl_m)+Zs*(vl_hat_m-vl_m))

  !call generate_fluctuations(Zp_out,Tn_p,Tn_hat_p,vn_p,vn_hat_p,FRn)
  !call generate_fluctuations(Zs_out,Tm_p,Tm_hat_p,vm_p,vm_hat_p,FRm)
  !call generate_fluctuations(Zs_out,Tl_p,Tl_hat_p,vl_p,vl_hat_p,FRl)

  FL_n = FLn/Zp
  FL_m = FLm/Zs
  FL_l = FLl/Zs

  !FR_n = FRn/Zp_out
  !FR_m = FRm/Zs_out
  !FR_l = FRl/Zs_out

  ! rotate to global
  call rotate_nml2xyz(n_n,n_m,n_l,FLn,FLm,FLl,FLx,FLy,FLz)
  !call rotate_nml2xyz(n_n,n_m,n_l,FRn,FRm,FRl,FRx,FRy,FRz)
  call rotate_nml2xyz(n_n,n_m,n_l,FL_n,FL_m,FL_l,FL_x,FL_y,FL_z)
  !call rotate_nml2xyz(n_n,n_m,n_l,FR_n,FR_m,FR_l,FR_x,FR_y,FR_z)

  fstar(1) = FLx
  fstar(2) = FLy
  fstar(3) = FLz
  fstar(4) = n_n(1)*FL_x
  fstar(5) = n_n(2)*FL_y
  fstar(6) = n_n(3)*FL_z
  fstar(7) = n_n(3)*FL_y+n_n(2)*FL_z
  fstar(8) = n_n(3)*FL_x+n_n(1)*FL_z
  fstar(9) = n_n(2)*FL_x+n_n(1)*FL_y

  !if (flux_method == 1 .and. &
  !    mesh%bctype(is,ie) == BC_IN .and. &
  !    mesh%fluxtype(is,ie)==1) then
  !  !fstar = (fR -fL)*0.5
  !  fstar(1) = dTx
  !  fstar(2) = dTy
  !  fstar(3) = dTz
  !  fstar(4) = n_n(1)*dVx
  !  fstar(5) = n_n(2)*dVy
  !  fstar(6) = n_n(3)*dVz
  !  fstar(7) = n_n(3)*dVy+n_n(2)*dVz
  !  fstar(8) = n_n(3)*dVx+n_n(1)*dVz
  !  fstar(9) = n_n(2)*dVx+n_n(1)*dVy
  !  fstar = 0.5*fstar
  !end if
100   fluxs(i+(is-1)*Nfp,:) = fstar

  end do
  end do

  do i = 1,Nvar
    fluxs(:,i) = mesh%Fscale(:,ie) * fluxs(:,i)
  end do

  !fluxes(((ie-1)*Nfp*Nfaces+1):(ie*Nfp*Nfaces),:) = fluxs

  !end do ! element

  !fstar = 0
end subroutine

subroutine time_weakening(y,z,y0,z0,rcrit,Vr,t0,cur_time,slip,Dc,mu_s,mu_d,mu_f)
  implicit none
  real(kind=RKIND) :: T,t0,f1,f2,dist,rcrit,Vr
  real(kind=RKIND) :: slip,Dc,mu_s,mu_d,mu_f
  real(kind=RKIND) :: y,z,y0,z0,cur_time

  dist = sqrt( (y-y0)**2 + (z-z0)**2 )
  if(dist<rcrit) then
    T=dist/Vr+0.081*rcrit/Vr*(1./(1.-(dist/rcrit)**2)-1.)
  else
    T = 1.0e9
  endif

  if(slip<Dc) then
    f1 = slip/Dc
  else
    f1 = 1.0
  endif
  if(cur_time<T) then
    f2 = 0.0
  elseif(cur_time<T+t0) then
    f2 = (cur_time-T)/t0
  else
    f2 = 1.0
  endif

  mu_f = mu_s-(mu_s-mu_d)*max(f1, f2)
end subroutine

! subroutine avg_face(mesh,u,unew,qi)
!   implicit none
!   type(meshvar) :: mesh
!   real(kind=RKIND),dimension(:,:) :: u,unew
!   real(kind=RKIND),dimension(:,:,:,:) :: qi
!   integer :: i,is,ie
!   real(kind=RKIND),dimension(Nvar) :: uL,uR,uLnew
!   integer :: neigh
!   integer :: direction
!   integer :: mpi_e,mpi_n!,n
!   integer :: flipped_index(Nfp)
!
!   do ie = 1,mesh%Nelem
!     do is = 1,Nfaces
!       do i = 1,Nfp
!         uL = u(mesh%vmapM(i,is,ie),:)
!
!         neigh = mesh%neigh(is,ie)
!         direction = mesh%direction(is,ie)
!         !if(neigh>0) uR=u(mesh%vmapP(i,is,ie),:)
!         if(neigh==ie) uR=0 ! absorbing
!         !if(neigh>0 .and. neigh .ne. ie) then
!         if(neigh>0) then
!           uR=u(mesh%vmapP(i,is,ie),:)
!         elseif (neigh==0) then
!           uR(:) = 0
!         elseif (neigh==-1) then
!           !print*,'mpi'
!           mpi_e = mesh%mpi_ibool(is,ie)
!           mpi_n = mesh%mpi_interface(4,is,ie)
!           !print*,'mpi_e=',mpi_e
!           !print*,'mpi_n=',mpi_n
!           flipped_index = mesh%flipped_index(:,direction)
!           uR(:) = qi(flipped_index(i),1:9,mpi_e,mpi_n)  ! mpi
!         else
!           uR(:) = 1e38
!         end if
!         if(neigh==ie) uR=0 ! absorbing
!
!         if (mesh%bctype(is,ie) == BC_FREE) then ! free surface
!           uR(1:3) = uL(1:3)
!           uR(4:9) = -uL(4:9)
!         end if
!
!         if (mesh%bctype(is,ie) == 0 .and. mesh%fluxtype(is,ie)==1) then
!           !fstar = (fR -fL)*0.5
!           unew(mesh%vmapM(i,is,ie),:) = uL !(uL+uR)*0.5
!         end if
!       end do
!     end do
!   end do
! end subroutine

! function norm3(a)
!   real(kind=RKIND) :: a(3)
!   real(kind=RKIND) :: norm3
!   norm3=dsqrt(a(1)**2+a(2)**2+a(3)**2)
! end
!
! subroutine Flux1(U,rho,cp,cs,F)
!   implicit none
!   real(kind=rkind),intent(in) :: U(:),rho,cp,cs
!   real(kind=rkind),intent(out) :: F(:)
!   real(kind=rkind) :: rrho,miu,lam,chi
!   rrho = 1d0/rho
!   miu = rho*cs**2
!   chi = rho*cp**2
!   lam = chi-2d0*miu
!
!   ! U=(rho*Vx,rho*Vy,rho*Vz,Exx,Eyy,Ezz,Eyz,Exz,Exy)
!   ! F=(Sxx,Sxy,Sxz,Vx,0,0,0,Vz,Vy)
!   F(1) = U(4) * chi + U(5) * lam + U(6) * lam
!   F(2) = U(9) * miu
!   F(3) = U(8) * miu
!   F(4) = U(1) * rrho
!   F(5) = 0
!   F(6) = 0
!   F(7) = 0
!   F(8) = U(3) * rrho
!   F(9) = U(2) * rrho
! end subroutine
!
! subroutine Flux2(U,rho,cp,cs,F)
!   implicit none
!   real(kind=rkind),intent(in) :: U(:),rho,cp,cs
!   real(kind=rkind),intent(out) :: F(:)
!   real(kind=rkind) :: rrho,miu,lam,chi
!   rrho = 1d0/rho
!   miu = rho*cs**2
!   chi = rho*cp**2
!   lam = chi-2d0*miu
!
!   ! U=(rho*Vx,rho*Vy,rho*Vz,Exx,Eyy,Ezz,Eyz,Exz,Exy)
!   ! F=(Sxy,Syy,Syz,0,Vy,0,Vz,0,Vx)
!   F(1) = U(9) * miu
!   F(2) = U(4) * lam + U(5) * chi + U(6) * lam
!   F(3) = U(7) * miu
!   F(4) = 0
!   F(5) = U(2) * rrho
!   F(6) = 0
!   F(7) = U(3) * rrho
!   F(8) = 0
!   F(9) = U(1) * rrho
! end subroutine
!
! subroutine Flux3(U,rho,cp,cs,F)
!   implicit none
!   real(kind=rkind),intent(in) :: U(:),rho,cp,cs
!   real(kind=rkind),intent(out) :: F(:)
!   real(kind=rkind) :: rrho,miu,lam,chi
!   rrho = 1d0/rho
!   miu = rho*cs**2
!   chi = rho*cp**2
!   lam = chi-2d0*miu
!
!   ! U=(rho*Vx,rho*Vy,rho*Vz,Exx,Eyy,Ezz,Eyz,Exz,Exy)
!   ! F=(Sxz,Syz,Szz,0,0,Vz,Vy,Vx,0)
!   F(1) = U(8) * miu
!   F(2) = U(7) * miu
!   F(3) = U(4) * lam + U(5) * lam + U(6) * chi
!   F(4) = 0
!   F(5) = 0
!   F(6) = U(3) * rrho
!   F(7) = U(2) * rrho
!   F(8) = U(1) * rrho
!   F(9) = 0
! end subroutine
! 
! subroutine strain2stress(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs,&
!                          sxx,syy,szz,syz,sxz,sxy)
!   implicit none
!   real(kind=RKIND),intent(in) :: exx,eyy,ezz,eyz,exz,exy,rho,cp,cs
!   real(kind=RKIND),intent(out) :: sxx,syy,szz,syz,sxz,sxy
!   !real(kind=rkind) :: rrho,miu,lam,chi
!   real(kind=rkind) :: miu,lam,chi
!   !rrho = 1d0/rho
!   miu = rho*cs**2
!   chi = rho*cp**2
!   lam = chi-2d0*miu
!
!   sxx = exx * chi + eyy * lam + ezz * lam
!   syy = exx * lam + eyy * chi + ezz * lam
!   szz = exx * lam + eyy * lam + ezz * chi
!   syz = eyz * miu
!   sxz = exz * miu
!   sxy = exy * miu
! end subroutine
!
! subroutine stress2strain(sxx,syy,szz,syz,sxz,sxy,rho,cp,cs,&
!                          exx,eyy,ezz,eyz,exz,exy)
!   implicit none
!   real(kind=RKIND),intent(in) :: sxx,syy,szz,syz,sxz,sxy,rho,cp,cs
!   real(kind=RKIND),intent(out) :: exx,eyy,ezz,eyz,exz,exy
!   !real(kind=rkind) :: rrho,miu,lam,chi
!   real(kind=rkind) :: miu,lam,chi,a,b
!   !rrho = 1d0/rho
!   miu = rho*cs**2
!   chi = rho*cp**2
!   lam = chi-2d0*miu
!   a = (chi+lam)/(chi**2+lam*chi-2.0*lam**2)
!   b = (-lam   )/(chi**2+lam*chi-2.0*lam**2)
!
!   exx = sxx * a + syy * b + szz * b
!   eyy = sxx * b + syy * a + szz * b
!   ezz = sxx * b + syy * b + szz * a
!   eyz = syz / miu
!   exz = sxz / miu
!   exy = sxy / miu
! end subroutine

! elemental subroutine prakash_cliff_fric(strength,sigma,V,Vstar,L,mu,dt)
!   implicit none
!   real(kind=RKIND), intent(inout):: strength !< shear strength (or friction)
!   real(kind=RKIND), intent(in)   :: V        !< slip velocity
!   real(kind=RKIND), intent(in)   :: sigma    !< normal traction
!   real(kind=RKIND), intent(in)   :: Vstar    !< referenz Velocity
!   real(kind=RKIND), intent(in)   :: L        !< referenz length
!   real(kind=RKIND), intent(in)   :: mu       !< friction coefficient
!   real(kind=RKIND), intent(in)   :: dt       !< time increment
!   real(kind=RKIND) :: expterm
!   expterm=exp(-(abs(V) + Vstar)*dt/L)
!   strength =   strength*expterm - max(0.,-mu*sigma)*(expterm-1.)
! end subroutine prakash_cliff_fric

!subroutine Return_Map(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs,depth,dt,&
!                      sxx,syy,szz,syz,sxz,sxy)
!@ subroutine Return_Map(sxx,syy,szz,syz,sxz,sxy,depth,dt)
!@   implicit none
!@   !real(kind=RKIND),intent(in) :: exx, eyy, ezz, exy, exz, eyz
!@   real(kind=RKIND),intent(inout) :: sxx, syy, szz, sxy, sxz, syz
!@   !real(kind=RKIND),intent(in) :: rho,cp,cs
!@   real(kind=RKIND),intent(in) :: depth,dt
!@
!@   real(kind=RKIND) :: sxx0, syy0, szz0, sxy0, sxz0, syz0, fluidpresh
!@   real(kind=RKIND) :: sm, sdxx, sdyy, sdzz, sdxy, sdxz, sdyz
!@   real(kind=RKIND) :: secinv, tau, taulim, decay, yldfac
!@   real(kind=RKIND) :: omeg
!@
!@   real(kind=RKIND) :: cohes, blkfric, angfric, tvisc, dist
!@   integer :: i, j, k
!@
!@   !real(kind=RKIND) :: lam,miu,chi
!@
!@   ! For tpv27
!@   cohes = 1.36e0
!@   blkfric = 0.1934
!@   tvisc = 0.03
!@
!@   !!! For tpv29
!@   !cohes = 1.18e6
!@   !blkfric = 0.1680
!@   !tvisc = 0.05
!@
!@   angfric = atan(blkfric)
!@
!@   !do k = nk1, nk2
!@   !  do j = nj1, nj2
!@   !    do i = ni1, ni2
!@         !sxx = Txx(i,j,k); syy = Tyy(i,j,k); szz = Tzz(i,j,k)
!@         !sxy = Txy(i,j,k); sxz = Txz(i,j,k); syz = Tyz(i,j,k)
!@         !depth = -z(i,j,k)
!@         !depth = -z
!@
!@         ! For tpv27
!@         fluidpresh = 1.0*9.8*depth ! rho*g*h
!@         if(depth<=15.0)then
!@           omeg = 1.0
!@         elseif(depth<=20.0)then
!@           omeg = (20.0 - depth)/5.0
!@         else
!@           omeg = 0.0
!@         endif
!@         !szz0 = min(-2670.0*9.8*steph/3.0, -2670.0*9.8*depth)
!@         !szz0 = -2670.0*9.8*depth
!@         szz0 = -2.670*9.8*depth ! in MPa
!@         syy0 = omeg*(0.926793*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
!@         sxx0 = omeg*(1.073206*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
!@         sxy0 = omeg*(0.169029*(szz0+fluidpresh))
!@         sxz0 = 0.0
!@         syz0 = 0.0
!@
!@         ! For tpv30
!@         !fluidpresh = depth*9.8e3
!@         !if(depth<=17000.0)then
!@         !  omeg = 1.0
!@         !elseif(depth<=22000.0)then
!@         !  omeg = (22000.0 - depth)/5000.0
!@         !else
!@         !  omeg = 0.0
!@         !endif
!@         !szz0 = min(-2670.0*9.8*steph/10.0, -2670.0*9.8*depth)
!@         !syy0 = omeg*(1.025837*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
!@         !sxx0 = omeg*(0.974162*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
!@         !sxy0 = omeg*(0.158649*(szz0+fluidpresh))
!@         !sxz0 = 0.0
!@         !syz0 = 0.0
!@
!@         !miu = rho*cs*cs
!@         !chi = rho*cp*cp
!@         !lam = chi-2d0*miu
!@
!@         !sxx = sxx + chi*exx + lam*eyy + lam*ezz
!@         !syy = syy + lam*exx + chi*eyy + lam*ezz
!@         !szz = szz + lam*exx + lam*eyy + chi*ezz
!@         !sxy = sxy + miu*exy
!@         !sxz = sxz + miu*exz
!@         !syz = syz + miu*eyz
!@
!@         sxx = sxx0 + sxx
!@         syy = syy0 + syy
!@         szz = szz0 + szz
!@         sxy = sxy0 + sxy
!@         sxz = sxz0 + sxz
!@         syz = syz0 + syz
!@
!@         sm = (sxx+syy+szz)/3.0
!@         sdxx = sxx - sm
!@         sdyy = syy - sm
!@         sdzz = szz - sm
!@         sdxy = sxy
!@         sdxz = sxz
!@         sdyz = syz
!@
!@         secinv = 0.5*(sdxx**2+sdyy**2+sdzz**2)+sdxy**2+sdxz**2+sdyz**2
!@         tau = sqrt(secinv)
!@         taulim = cohes*cos(angfric) - (sm+fluidpresh)*sin(angfric)
!@         taulim = max(0.0, taulim)
!@
!@         if(tau .gt. taulim) then
!@           !decay = exp(-stept/tvisc)
!@           decay = exp(-dt/tvisc)
!@           yldfac = decay + (1.0-decay)*taulim/tau
!@           sxx = sdxx*yldfac + sm
!@           syy = sdyy*yldfac + sm
!@           szz = sdzz*yldfac + sm
!@           sxy = sdxy*yldfac
!@           sxz = sdxz*yldfac
!@           syz = sdyz*yldfac
!@         endif
!@ 
!@         sxx = sxx - sxx0
!@         syy = syy - syy0
!@         szz = szz - szz0
!@         sxy = sxy - sxy0
!@         sxz = sxz - sxz0
!@         syz = syz - syz0
!@         !Txx(i,j,k) = sxx; Tyy(i,j,k)= syy;  Tzz(i,j,k) = szz
!@         !Txy(i,j,k) = sxy; Txz(i,j,k)= sxz;  Tyz(i,j,k) = syz
!@   !    enddo
!@   !  enddo
!@   !enddo
!@ endsubroutine
!@
!@ subroutine update_plastic(mesh,u)
!@   implicit none
!@   type(meshvar) :: mesh
!@   real(kind=RKIND),dimension(:,:) :: u ! (Np*Nelem,Nvar)
!@   real(kind=RKIND) :: dt,depth,rho,cp,cs
!@   real(kind=RKIND) :: exx,eyy,ezz,eyz,exz,exy
!@   real(kind=RKIND) :: sxx,syy,szz,syz,sxz,sxy
!@   integer :: ie,i
!@
!@   dt = mesh%deltat
!@
!@   do ie = 1,mesh%Nelem
!@     !do j = 1,Np
!@     !  iv(j) = j+(ie-1)*Np
!@     !end do
!@     cp = mesh%vp(ie)
!@     cs = mesh%vs(ie)
!@     rho = mesh%rho(ie)
!@     !rrho = 1d0/rho
!@
!@     do i = 1,Np
!@
!@       depth = -mesh%vz(i+(ie-1)*Np)
!@       exx = u(i+(ie-1)*Np,4)
!@       eyy = u(i+(ie-1)*Np,5)
!@       ezz = u(i+(ie-1)*Np,6)
!@       eyz = u(i+(ie-1)*Np,7)
!@       exz = u(i+(ie-1)*Np,8)
!@       exy = u(i+(ie-1)*Np,9)
!@
!@       call strain2stress(exx,eyy,ezz,eyz,exz,exy,rho,cp,cs, &
!@                          sxx,syy,szz,syz,sxz,sxy)
!@
!@       call Return_Map(sxx,syy,szz,syz,sxz,sxy,depth,dt)
!@
!@       call stress2strain(sxx,syy,szz,syz,sxz,sxy,rho,cp,cs, &
!@                          exx,eyy,ezz,eyz,exz,exy)
!@
!@       u(i+(ie-1)*Np,4) = exx
!@       u(i+(ie-1)*Np,5) = eyy
!@       u(i+(ie-1)*Np,6) = ezz
!@       u(i+(ie-1)*Np,7) = eyz
!@       u(i+(ie-1)*Np,8) = exz
!@       u(i+(ie-1)*Np,9) = exy
!@     end do
!@   end do ! element
!@
!@ end subroutine

subroutine init_wave(mesh,u)
  implicit none
  type(meshvar) :: mesh
  real(kind=RKIND),dimension(:,:) :: u ! (Np*Nelem,Nvar)
  real(kind=RKIND) :: srcx,srcy,srcz,a
  real(kind=RKIND) :: xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=RKIND),allocatable,dimension(:) :: r,amp
  integer :: Nelem

  if (initial_condition_wave .ne. 1) return

  Nelem = mesh%Nelem

  allocate(r(Np*Nelem))
  allocate(amp(Np*Nelem))

  xmin = minval(mesh%vx)
  xmax = maxval(mesh%vx)
  ymin = minval(mesh%vy)
  ymax = maxval(mesh%vy)
  zmin = minval(mesh%vz)
  zmax = maxval(mesh%vz)

  !print*, 'x=(',sngl(xmin),',',sngl(xmax),')'
  !print*, 'y=(',sngl(ymin),',',sngl(ymax),')'
  !print*, 'z=(',sngl(zmin),',',sngl(zmax),')'

  srcx = (xmin+xmax)/2.0
  srcy = (ymin+ymax)/2.0
  srcz = (zmin+zmax)/2.0
  srcx = src_loc(1)
  srcy = src_loc(2)
  srcz = src_loc(3)

  r = (mesh%vx-srcx)**2+(mesh%vy-srcy)**2+(mesh%vz-srcz)**2
  a = src_gaussian_width
  amp = src_m0*exp(-r/a**2)

  ! rhoVx rhoVy rhoVz Exx Eyy Ezz Eyz Exz Exy
  u(:,1) = 0
  u(:,2) = 0
  u(:,3) = 0
  u(:,4) = amp*src_mxx
  u(:,5) = amp*src_myy
  u(:,6) = amp*src_mzz
  u(:,7) = amp*src_myz
  u(:,8) = amp*src_mxz
  u(:,9) = amp*src_mxy

  deallocate(r,amp)

end subroutine

!subroutine write_wave(u,rank,flag)
!  implicit none
!  real(kind=rkind) :: u(:,:)
!  integer,intent(in) :: rank,flag
!  character(len=128) :: filename
!  !write(filename,'(a,i6.6,a)') 'data/wave_it',it,'.bin'
!  !open(100,file=trim(filename),access='direct',form='unformatted',recl=sizeof(1.0)*size(u),status='replace')
!  !write(100,rec=1) sngl(u)
!
!  write(filename,'(a,i6.6)') 'data/wave',rank
!  if (flag == 0) then
!      open(100,file=trim(filename),access='stream',form='unformatted',status='replace',action='write')
!  else
!      open(100,file=trim(filename),access='stream',form='unformatted',status='unknown',position='append')
!  end if
!  write(100) sngl(u)
!  close(100)
!
!end subroutine

function Gt_func(t,rise_time) result (g)
  real(kind=RKIND) :: g
  real(kind=RKIND) :: t,rise_time
  real(kind=RKIND) :: t1,t2
  if (t < rise_time) then
    t1 = t-rise_time
    t2 = t-2.0*rise_time
    g = exp(t1**2/t/t2)
  else
    g = 1.0
  end if

end function

subroutine write_wave_vtk(mesh,u)
  implicit none
  type(meshvar) :: mesh
  integer :: myrank
  integer :: ie,i,j,k
  character(len=256) :: filename
  real(kind=rkind),dimension(:,:) :: u
  real(kind=rkind),dimension(:),allocatable :: tmp

  myrank = mesh%rank

  i = export_wave_component
  allocate(tmp(mesh%Nelem))
  do ie = 1,mesh%Nelem
    j=(ie-1)*Np+1
    k=(ie-1)*Np+Np
    tmp(ie) = sum(u(j:k,i))/dble(Np)
    if(i<=3) then
      tmp(ie) = tmp(ie)/mesh%rho(ie)
    end if
  end do
  write(filename,'(a,a,i6.6,a)') trim(mesh_dir),'/wave',myrank,'.vtk'
  write(*,*) "write ", trim(filename)
  call writeVtkTetraMeshRealdata(filename, mesh%elem, mesh%coord, tmp)
  deallocate(tmp)
end subroutine

end module
