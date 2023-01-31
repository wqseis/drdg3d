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

module mod_init_fault

  use mod_para,   only : RKIND,           &
                         Nfp, Np, Nfaces, &
                         BC_FAULT,        &
                         friction_law,    &
                         thermalpressure, &
                         mesh_dir
  use mod_types,  only : meshvar
  use mod_rotate, only : rotate_xyz2nml
  use netcdf

  implicit none

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

! subroutine fault_init_tpv5(mesh)
!   implicit none
!   type(meshvar) :: mesh
!   real(rkind) :: x,y,z,xc,yc,zc,asp_size
!   integer :: is,ie,ief,i,j!,k
!   integer :: myrank
!   real(rkind) :: sxx,syy,szz,sxy,sxz,syz
!   real(rkind) :: vec_n(3),vec_m(3),vec_l(3),Tx,Ty,Tz,Tn,Tm,Tl
!   !real(rkind) :: damp,r,r1,r2
!   !real(rkind) :: Pf,C0,depth
!   real(rkind) :: depth
!   !real(rkind) :: xface(Nfp),yface(Nfp)
!   !character(len=128) :: filename
!   !real*8,parameter :: b22 =  0.926793
!   !real*8,parameter :: b33 =  1.073206
!   !real*8,parameter :: b23 = -0.169029
!   !real(rkind) :: T_global(3,3),T_local(3,3)
!   !real(rkind) :: mat_L2G(3,3),e_local(3,3),e_global(3,3)
!   !real(rkind) :: h1,h2,Dc
!   !logical :: isfault
! 
!   myrank = mesh%rank
!   ! ------------------------------------------------------------
!   ! stress and friction initialization 
! 
!   do ief = 1,mesh%nfault_elem
!     ie = mesh%fault2wave(ief)
!     do is = 1,Nfaces
!       if(mesh%bctype(is,ie) >= BC_FAULT) then
!         xc = sum(mesh%vx(mesh%vmapM(:,is,ie)))/dble(Nfp)
!         yc = sum(mesh%vy(mesh%vmapM(:,is,ie)))/dble(Nfp)
!         zc = sum(mesh%vz(mesh%vmapM(:,is,ie)))/dble(Nfp)
!         do i = 1,Nfp
!           j = i+(is-1)*Nfp
!           x = mesh%vx(mesh%vmapM(i,is,ie))
!           y = mesh%vy(mesh%vmapM(i,is,ie))
!           z = mesh%vz(mesh%vmapM(i,is,ie))
! 
!           depth = abs(z)
! 
!           vec_n = (/mesh%nx(j,ie),mesh%ny(j,ie),mesh%nz(j,ie)/)
!           vec_m = (/mesh%mx(j,ie),mesh%my(j,ie),mesh%mz(j,ie)/)
!           vec_l = (/mesh%lx(j,ie),mesh%ly(j,ie),mesh%lz(j,ie)/)
! 
!           !if(masternode) print*,'set tpv5/6 frictions'
! 
!           mesh%mu_s (i,is,ief) = 0.677
!           mesh%mu_d (i,is,ief) = 0.525
!           mesh%Dc   (i,is,ief) = 0.4
!           mesh%C0   (i,is,ief) = 0.0e0
! 
!           sxx = -120e0
!           syy = 0
!           szz = 0
!           sxy = -70e0
!           sxz = 0
!           syz = 0
! 
!           asp_size = 1.5e0
!           if ( abs(y-0e0)<=asp_size .and. abs(z+7.5e0)<=asp_size ) then
!             sxy = -81.6e0
!           end if
! 
!           ! rectangle aspersity
!           asp_size = 1.5e0
!           if ( abs(y-0e0)<=asp_size .and. abs(z+7.5e0)<=asp_size ) then
!             sxy = -81.6e0
!           end if
!           ! lower asperity in the right
!           if ( abs(y-7.5e0)<=asp_size .and. abs(z+7.5e0)<=asp_size ) then
!             sxy = -62e0
!           end if
!           ! higher asperity in the left
!           if ( abs(y+7.5e0)<=asp_size .and. abs(z+7.5e0)<=asp_size ) then
!             sxy = -78e0
!           end if
! 
!           Tx = sxx*vec_n(1)+sxy*vec_n(2)+sxz*vec_n(3)
!           Ty = sxy*vec_n(1)+syy*vec_n(2)+syz*vec_n(3)
!           Tz = sxz*vec_n(1)+syz*vec_n(2)+szz*vec_n(3)
!           call rotate_xyz2nml(vec_n,vec_m,vec_l,Tx,Ty,Tz,Tn,Tm,Tl)
! 
!           mesh%tau0n(i,is,ief) = Tn
!           mesh%tau0m(i,is,ief) = Tm
!           mesh%tau0l(i,is,ief) = Tl
!         enddo
!       else
!         mesh%mu_s(:,is,ief) = 1e4
!         mesh%C0  (:,is,ief) = 1e9
!       endif
!     enddo
!   enddo
! 
!   mesh%stress = mesh%tau0n
! 
! end subroutine

subroutine fault_init_external(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: is,ie,ief,i,j
  integer :: myrank
  integer :: ncid,varid,ierr
  character(len=128) :: filename
  real(kind=rkind) :: vec_n(3),vec_m(3),vec_l(3)
  real(kind=rkind) :: Tx,Ty,Tz,Tn,Tm,Tl
  real(kind=rkind) :: dTx,dTy,dTz,dTn,dTm,dTl
  real(kind=rkind),allocatable,dimension(:,:,:) :: Tx0,Ty0,Tz0
  real(kind=rkind),allocatable,dimension(:,:,:) :: mu_s,mu_d,Dc,C0
  real(kind=rkind),allocatable,dimension(:,:,:) :: dTx0,dTy0,dTz0
  real(kind=rkind) :: c(3),cp,p(3),v1(3),v2(3),v3(3)
  ! rate state
  real(kind=rkind),allocatable,dimension(:,:,:) :: a,b,Vw,state
  ! thermal pressurization
  real(kind=rkind),allocatable,dimension(:,:,:) :: TP_hy

  integer :: FToV(4,3)

  ! point outward
  ! Face 1: 1,3,2
  ! Face 2: 1,2,4
  ! Face 3: 2,3,4
  ! Face 4: 1,4,3
  FToV = reshape((/1,1,2,1,3,2,3,4,2,4,4,3/),shape(FToV))

  myrank = mesh%rank

  !allocate(temp(3,Nfaces,mesh%nfault_elem))
  allocate(Tx0 (3,Nfaces,mesh%nfault_elem))
  allocate(Ty0 (3,Nfaces,mesh%nfault_elem))
  allocate(Tz0 (3,Nfaces,mesh%nfault_elem))
  allocate(dTx0(3,Nfaces,mesh%nfault_elem))
  allocate(dTy0(3,Nfaces,mesh%nfault_elem))
  allocate(dTz0(3,Nfaces,mesh%nfault_elem))
  allocate(mu_s(3,Nfaces,mesh%nfault_elem))
  allocate(mu_d(3,Nfaces,mesh%nfault_elem))
  allocate(Dc  (3,Nfaces,mesh%nfault_elem))
  allocate(C0  (3,Nfaces,mesh%nfault_elem))
  ! rate state
  allocate(a    (3,Nfaces,mesh%nfault_elem))
  allocate(b    (3,Nfaces,mesh%nfault_elem))
  allocate(Vw   (3,Nfaces,mesh%nfault_elem))
  allocate(state(3,Nfaces,mesh%nfault_elem))
  ! thermal pressurization
  allocate(TP_hy(3,Nfaces,mesh%nfault_elem))

  write(filename,'(a,a,i6.6,a)') trim(mesh_dir),'/meshVar',myrank,'.nc'
  ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid)

  ierr = nf90_inq_varid(ncid,'Tx0',varid)
  call check2(ierr,'inq_varid Tx0')
  ierr = nf90_get_var(ncid,varid,Tx0)
  call check2(ierr,'get_var Tx0')

  ierr = nf90_inq_varid(ncid,'Ty0',varid)
  call check2(ierr,'inq_varid Ty0')
  ierr = nf90_get_var(ncid,varid,Ty0)
  call check2(ierr,'get_var Ty0')

  ierr = nf90_inq_varid(ncid,'Tz0',varid)
  call check2(ierr,'inq_varid Tz0')
  ierr = nf90_get_var(ncid,varid,Tz0)
  call check2(ierr,'get_var Tz0')

  ierr = nf90_inq_varid(ncid,'dTx0',varid)
  call check2(ierr,'inq_varid dTx0')
  ierr = nf90_get_var(ncid,varid,dTx0)
  call check2(ierr,'get_var dTx0')

  ierr = nf90_inq_varid(ncid,'dTy0',varid)
  call check2(ierr,'inq_varid dTy0')
  ierr = nf90_get_var(ncid,varid,dTy0)
  call check2(ierr,'get_var dTy0')

  ierr = nf90_inq_varid(ncid,'dTz0',varid)
  call check2(ierr,'inq_varid dTz0')
  ierr = nf90_get_var(ncid,varid,dTz0)
  call check2(ierr,'get_var dTz0')

  ierr = nf90_inq_varid(ncid,'mu_s',varid)
  call check2(ierr,'inq_varid mu_s')
  ierr = nf90_get_var(ncid,varid,mu_s)
  call check2(ierr,'get_var mu_s')

  ierr = nf90_inq_varid(ncid,'mu_d',varid)
  call check2(ierr,'inq_varid mu_d')
  ierr = nf90_get_var(ncid,varid,mu_d)
  call check2(ierr,'get_var mu_d')

  ierr = nf90_inq_varid(ncid,'Dc',varid)
  call check2(ierr,'inq_varid Dc')
  ierr = nf90_get_var(ncid,varid,Dc)
  call check2(ierr,'get_var Dc')

  ierr = nf90_inq_varid(ncid,'C0',varid)
  call check2(ierr,'inq_varid C0')
  ierr = nf90_get_var(ncid,varid,C0)
  call check2(ierr,'get_var C0')

  if (friction_law == 1 .or. &
      friction_law == 2 .or. &
      friction_law == 3) then
    ! rate state
    ierr = nf90_inq_varid(ncid,'a',varid)
    call check2(ierr,'inq_varid a')
    ierr = nf90_get_var(ncid,varid,a)
    call check2(ierr,'get_var a')

    ierr = nf90_inq_varid(ncid,'b',varid)
    call check2(ierr,'inq_varid b')
    ierr = nf90_get_var(ncid,varid,b)
    call check2(ierr,'get_var b')

    ierr = nf90_inq_varid(ncid,'state',varid)
    call check2(ierr,'inq_varid state')
    ierr = nf90_get_var(ncid,varid,state)
    call check2(ierr,'get_var state')

    if (friction_law == 3) then
      ierr = nf90_inq_varid(ncid,'Vw',varid)
      call check2(ierr,'inq_varid Vw')
      ierr = nf90_get_var(ncid,varid,Vw)
      call check2(ierr,'get_var Vw')
    end if

  end if

  if (thermalpressure == 1) then
    ierr = nf90_inq_varid(ncid,'TP_hy',varid)
    call check2(ierr,'inq_varid TP_hy')
    ierr = nf90_get_var(ncid,varid,TP_hy)
    call check2(ierr,'get_var TP_hy')
  end if

  do ief = 1,mesh%nfault_elem
    ie = mesh%fault2wave(ief)
    do is = 1,Nfaces
      ! naive triangular interpolation,
      ! will change to barycentric method later
      if (mesh%bctype(is,ie) >= BC_FAULT) then

        v1=mesh%coord(1:3,mesh%elem(FtoV(is,1),ie))
        v2=mesh%coord(1:3,mesh%elem(FtoV(is,2),ie))
        v3=mesh%coord(1:3,mesh%elem(FtoV(is,3),ie))

        do i = 1,Nfp
          p(1)=mesh%vx(mesh%vmapM(i,is,ie))
          p(2)=mesh%vy(mesh%vmapM(i,is,ie))
          p(3)=mesh%vz(mesh%vmapM(i,is,ie))

          c=mu_s(1:3,is,ief)
          cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
          mesh%mu_s(i,is,ief)=cp

          c=mu_d(1:3,is,ief)
          cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
          mesh%mu_d(i,is,ief)=cp

          c=Dc(1:3,is,ief)
          cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
          mesh%Dc(i,is,ief)=cp

          c=C0(1:3,is,ief)
          cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
          mesh%C0(i,is,ief)=cp

          c=Tx0(1:3,is,ief)
          Tx=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)

          c=Ty0(1:3,is,ief)
          Ty=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)

          c=Tz0(1:3,is,ief)
          Tz=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)

          c=dTx0(1:3,is,ief)
          dTx=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)

          c=dTy0(1:3,is,ief)
          dTy=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)

          c=dTz0(1:3,is,ief)
          dTz=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)

          if (friction_law == 1 .or. &
              friction_law == 2 .or. &
              friction_law == 3) then
            ! rate state
            c=a(1:3,is,ief)
            cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
            mesh%a(i,is,ief)=cp

            c=b(1:3,is,ief)
            cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
            mesh%b(i,is,ief)=cp

            c=state(1:3,is,ief)
            cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
            mesh%state(i,is,ief)=cp

            if (friction_law == 3) then
              c=Vw(1:3,is,ief)
              cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
              mesh%Vw(i,is,ief)=cp
            end if

            if (thermalpressure == 1) then
              c=TP_hy(1:3,is,ief)
              cp=tri_interp_dist(v1,v2,v3,c(1),c(2),c(3),p)
              mesh%TP_hy(i,is,ief)=cp
            end if

          end if

          j = i+(is-1)*Nfp
          vec_n = (/mesh%nx(j,ie),mesh%ny(j,ie),mesh%nz(j,ie)/)
          vec_m = (/mesh%mx(j,ie),mesh%my(j,ie),mesh%mz(j,ie)/)
          vec_l = (/mesh%lx(j,ie),mesh%ly(j,ie),mesh%lz(j,ie)/)


          call rotate_xyz2nml(vec_n,vec_m,vec_l,Tx,Ty,Tz,Tn,Tm,Tl)
          call rotate_xyz2nml(vec_n,vec_m,vec_l,dTx,dTy,dTz,dTn,dTm,dTl)

          mesh%tau0n(i,is,ief) = Tn
          mesh%tau0m(i,is,ief) = Tm
          mesh%tau0l(i,is,ief) = Tl
          mesh%dtau0n(i,is,ief) = dTn
          mesh%dtau0m(i,is,ief) = dTm
          mesh%dtau0l(i,is,ief) = dTl

          ! write initial stress
          mesh%stress1(i,is,ief) = Tm+0*dTm
          mesh%stress2(i,is,ief) = Tl+0*dTl
          mesh%stress(i,is,ief) = sqrt((Tm+0*dTm)**2+(Tl+0*dTl)**2)
          mesh%sigma (i,is,ief) = Tn+0*dTn
        end do

      else! not fault
        mesh%mu_s(:,is,ief) = 1e4
        mesh%C0  (:,is,ief) = 1e9
      end if
    end do
  end do

  ierr = nf90_close(ncid)

  deallocate(Tx0,Ty0,Tz0)
  deallocate(mu_s,mu_d,Dc,C0)
  ! rate state
  deallocate(a,b,Vw,state)
  ! thermal pressurization
  deallocate(TP_hy)
end subroutine

function tri_interp_dist(v1,v2,v3,c1,c2,c3,p) result (cp)
  implicit none

  real(kind=rkind),dimension(3) :: v1,v2,v3
  real(kind=rkind) :: c1,c2,c3,cp
  real(kind=rkind),dimension(3) :: p
  real(kind=rkind) :: r1,r2,r3
  real(kind=rkind) :: w1,w2,w3

  r1=sqrt((p(1)-v1(1))**2+(p(2)-v1(2))**2+(p(3)-v1(3))**2)
  r2=sqrt((p(1)-v2(1))**2+(p(2)-v2(2))**2+(p(3)-v2(3))**2)
  r3=sqrt((p(1)-v3(1))**2+(p(2)-v3(2))**2+(p(3)-v3(3))**2)

  w1=1d0/(r1+1d-300)
  w2=1d0/(r2+1d-300)
  w3=1d0/(r3+1d-300)

  cp=(w1*c1+w2*c2+w3*c3)/(w1+w2+w3)

end function

end module
