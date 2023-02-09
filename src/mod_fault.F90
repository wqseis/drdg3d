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

module mod_fault
  use mod_para,       only : RKIND,             &
                             Np, Nfp, Nfaces,   &
                             problem,           &
                             BC_FREE, BC_FAULT
  use mod_types,      only : meshvar
  use mod_init_fault, only : fault_init_external, &
                             fault_init_tpv3
  !use mod_rotate,     only : rotate_xyz2nml

  implicit none

contains

subroutine fault_init(mesh)
  implicit none
  type(meshvar) :: mesh
  !real(rkind) :: x,y,z,xc,yc,zc!,asp_size
  integer :: ie,is,k
  !integer :: nfault_elem,myrank
  !real(rkind) :: sxx,syy,szz,sxy,sxz,syz
  !real(rkind) :: vec_n(3),vec_m(3),vec_l(3),Tx,Ty,Tz,Tn,Tm,Tl
  !real(rkind) :: damp,r,r1,r2
  !real(rkind) :: Pf,C0,omeg,depth
  !real(rkind) :: xface(Nfp),yface(Nfp)
  !character(len=128) :: filename
  !real*8,parameter :: b22 =  0.926793
  !real*8,parameter :: b33 =  1.073206
  !real*8,parameter :: b23 = -0.169029
  !real(rkind) :: T_global(3,3),T_local(3,3)
  !real(rkind) :: mat_L2G(3,3),e_local(3,3),e_global(3,3)
  !real(rkind) :: h1,h2,angle,Dc,Dc0,p_ratio,sigma_v
  !integer :: m,n
  logical :: isfault

  !myrank = mesh%rank

  ! ------------------------------------------------------------
  ! fault variables initialization 

  mesh%nfault_elem = 0
  do ie = 1,mesh%nelem
  do is = 1,Nfaces
  if (mesh%bctype(is,ie) >= BC_FAULT) then
      mesh%nfault_elem = mesh%nfault_elem + 1
  end if
  end do
  end do

  mesh%nfault_elem = 0
  do ie = 1,mesh%nelem
    isfault=.false.
    do is = 1,Nfaces
      if (mesh%bctype(is,ie) >= BC_FAULT) then
        isfault=.true.
      end if
    end do
    if (isfault) then
        mesh%nfault_elem = mesh%nfault_elem + 1
    end if
  end do

  mesh%nfault_face = 0
  do ie = 1,mesh%nelem
    do is = 1,Nfaces
      if (mesh%bctype(is,ie) >= BC_FAULT) then
        mesh%nfault_face = mesh%nfault_face + 1
      end if
    end do
  end do

  ! count free faces for ground surface output
  mesh%nfree_face = 0
  do ie = 1,mesh%nelem
    do is = 1,Nfaces
      if (mesh%bctype(is,ie) == BC_FREE) then
        mesh%nfree_face = mesh%nfree_face + 1
      end if
    end do
  end do

  !print*,'rank=',myrank,'nfault=',mesh%nfault_elem,&
  !'nfree=',mesh%nfree_face

  allocate(mesh%fault2wave(mesh%nfault_elem))
  allocate(mesh%wave2fault(mesh%nelem))

  mesh%wave2fault(:) = 0

  k = 0
  do ie = 1,mesh%nelem
    isfault=.false.
    do is = 1,Nfaces
      if (mesh%bctype(is,ie) >= BC_FAULT) then
        isfault=.true.
      end if
    end do
    if (isfault) then
      k = k + 1
      mesh%fault2wave(k) = ie
      mesh%wave2fault(ie) = k
    end if
  end do

  !allocate(mesh%tau0n   (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%tau0m   (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%tau0l   (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%stress  (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%sigma   (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%Slip    (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%mSlip   (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%tSlip   (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%sliprate(Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%sliprate1(Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%sliprate2(Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%mu_s    (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%mu_d    (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%Dc      (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%C0      (Nfp,Nfaces,mesh%Nelem))
  !allocate(mesh%ruptime (Nfp,Nfaces,mesh%Nelem))

  allocate(mesh%tau0n   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%tau0m   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%tau0l   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%dtau0n   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%dtau0m   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%dtau0l   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%stress  (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%stress1 (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%stress2 (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%sigma   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%Slip    (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%Slip1   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%Slip2   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%mSlip1  (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%tSlip1  (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%mSlip2  (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%tSlip2  (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%sliprate(Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%sliprate1(Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%sliprate2(Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%mu_s    (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%mu_d    (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%Dc      (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%C0      (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%ruptime (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%peakrate(Nfp,Nfaces,mesh%nfault_elem))

  ! rate state
  allocate(mesh%a    (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%b    (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%Vw   (Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%state(Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%hstate(Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%tstate(Nfp,Nfaces,mesh%nfault_elem))
  allocate(mesh%mstate(Nfp,Nfaces,mesh%nfault_elem))
  ! thermal pressurization
  allocate(mesh%TP_hy(Nfp,Nfaces,mesh%nfault_elem))


  mesh%hstate(:,:,:) = 0
  mesh%tstate(:,:,:) = 0
  mesh%mstate(:,:,:) = 0

  mesh%slip    (:,:,:) = 0
  mesh%slip1   (:,:,:) = 0
  mesh%slip2   (:,:,:) = 0
  mesh%mslip1   (:,:,:) = 0
  mesh%tslip1   (:,:,:) = 0
  mesh%mslip2   (:,:,:) = 0
  mesh%tslip2   (:,:,:) = 0
  mesh%sliprate(:,:,:) = 0
  mesh%sliprate1(:,:,:) = 0
  mesh%sliprate2(:,:,:) = 0

  mesh%ruptime(:,:,:) = -1
  mesh%peakrate(:,:,:) = 0


  !call fault_init_tpv5(mesh)
  call fault_init_external(mesh)

  if (&
      trim(adjustl(problem)) .eq. 'tpv3' .or. &
      trim(adjustl(problem)) .eq. 'TPV3' ) then
    call fault_init_tpv3(mesh)
  end if

  ! ------------------------------------------------------------

!@    nfault_elem = 0
!@    !do ie = 1,mesh%Nelem
!@    do ief = 1,mesh%nfault_elem
!@        ie = mesh%fault2wave(ief)
!@        do is = 1,Nfaces
!@            if(mesh%bctype(is,ie) >= BC_FAULT) then
!@                xc = sum(mesh%vx(mesh%vmapM(:,is,ie)))/dble(Nfp)
!@                yc = sum(mesh%vy(mesh%vmapM(:,is,ie)))/dble(Nfp)
!@                zc = sum(mesh%vz(mesh%vmapM(:,is,ie)))/dble(Nfp)
!@                do i = 1,Nfp
!@                    j = i+(is-1)*Nfp
!@                    x = mesh%vx(mesh%vmapM(i,is,ie))
!@                    y = mesh%vy(mesh%vmapM(i,is,ie))
!@                    z = mesh%vz(mesh%vmapM(i,is,ie))
!@
!@                    depth = abs(z)
!@
!@                    vec_n = (/mesh%nx(j,ie),mesh%ny(j,ie),mesh%nz(j,ie)/)
!@                    vec_m = (/mesh%mx(j,ie),mesh%my(j,ie),mesh%mz(j,ie)/)
!@                    vec_l = (/mesh%lx(j,ie),mesh%ly(j,ie),mesh%lz(j,ie)/)
!@
!@!!!#if defined (TPV5) || defined (TPV6) || defined (WENCHUAN)
!@                    if ( &
!@                        trim(adjustl(problem)) .eq. 'tpv5' .or. &
!@                        trim(adjustl(problem)) .eq. 'tpv6' .or. &
!@                        trim(adjustl(problem)) .eq. 'wenchuan' ) then
!@
!@                    !if(masternode) print*,'set tpv5/6 frictions'
!@
!@                    mesh%mu_s (i,is,ief) = 0.677
!@                    mesh%mu_d (i,is,ief) = 0.525
!@                    mesh%Dc   (i,is,ief) = 0.4
!@                    mesh%C0   (i,is,ief) = 0.0e0
!@                    !mesh%tau_n(i,is,ie) = -120e6
!@                    !mesh%tau_0(i,is,ie) = -70e6
!@
!@                    sxx = -120e0
!@                    syy = 0
!@                    szz = 0
!@                    sxy = -70e0
!@                    sxz = 0
!@                    syz = 0
!@
!@                    !%if (.false.) then
!@                    !%! circle asperity
!@                    !%r1 = 1.69e0
!@                    !%r2 = 3e0
!@                    !%asp_size = r2*(1.0+1e-15)
!@                    !%r = sqrt((y+0e0)**2+(z+7.5e0)**2)
!@                    !%if ( r < asp_size ) then
!@                    !%    if (r<r1) then
!@                    !%        damp = 1d0
!@                    !%    else
!@                    !%        damp = 0.5*(1+cos(pi*(r-r1)/(r2-r1)))
!@                    !%    end if
!@                    !%    sxy = -70e0 - damp * 11.6e0
!@                    !%    sxy = +70e0 + damp * 11.6e0
!@                    !%endif
!@                    !%end if
!@
!@                    asp_size = 1.5e0
!@                    if ( abs(y-0e0)<=asp_size .and. abs(z+7.5e0)<=asp_size ) then
!@                        sxy = -81.6e0
!@                    end if
!@
!@                    end if
!@!!#endif
!@
!@!!!#ifdef TPV5
!@                    if ( trim(adjustl(problem)) .eq. 'tpv5' ) then
!@                    ! rectangle aspersity
!@                    !if (.true.) then
!@                    asp_size = 1.5e0
!@                    if ( abs(y-0e0)<=asp_size .and. abs(z+7.5e0)<=asp_size ) then
!@                        sxy = -81.6e0
!@                    end if
!@                    !if (.false.) then
!@                    ! lower asperity in the right
!@                    if ( abs(y-7.5e0)<=asp_size .and. abs(z+7.5e0)<=asp_size ) then
!@                        sxy = -62e0
!@                    end if
!@                    ! higher asperity in the left
!@                    if ( abs(y+7.5e0)<=asp_size .and. abs(z+7.5e0)<=asp_size ) then
!@                        sxy = -78e0
!@                    end if
!@                    !end if
!@                    !end if
!@                    end if
!@!!!#endif
!@
!@!!!#ifdef TPV24
!@                    if ( trim(adjustl(problem)) .eq. 'tpv24' ) then
!@                    Pf = 1.0*9.8*(-z) ! In MPa
!@                    szz = 2.670*9.8*(z)
!@                    if (z>-15.6) then
!@                        syy = b22*(szz+Pf)-Pf
!@                        sxx = b33*(szz+Pf)-Pf
!@                        sxy = b23*(szz+Pf)
!@                    else
!@                        syy = szz
!@                        sxx = szz
!@                        sxy = 0
!@                    end if
!@                    sxy = -sxy
!@                    syz = 0
!@                    sxz = 0
!@
!@                    if (z>-4.0) then
!@                        C0 = 0.3+0.675*(4.0+z)
!@                    else
!@                        C0 = 0.3
!@                    end if
!@
!@                    mesh%mu_s (i,is,ief) = 0.18
!@                    mesh%mu_d (i,is,ief) = 0.12
!@                    mesh%Dc   (i,is,ief) = 0.30
!@                    mesh%C0   (i,is,ief) = C0
!@                    end if
!@!!!#endif
!@
!@!!!#ifdef TPV26
!@                    if ( trim(adjustl(problem)) .eq. 'tpv26' ) then
!@                    if(z>-15.0)then
!@                        omeg = 1.0
!@                    elseif(z>-20.0)then
!@                        omeg = (20.0+z)/5.0
!@                    else
!@                        omeg = 0.0
!@                    endif
!@                    Pf = 1.0*9.8*(-z) ! In MPa
!@                    szz = 2.670*9.8*(z)
!@                    syy = omeg*(0.926793*(szz+Pf)-Pf) + (1.0-omeg)*szz !\sigma_11
!@                    sxx = omeg*(1.073206*(szz+Pf)-Pf) + (1.0-omeg)*szz !\sigma_33
!@                    sxy =-omeg*(0.169029*(szz+Pf))                     !\sigma_13
!@                    sxy = -sxy
!@                    syz = 0
!@                    sxz = 0
!@
!@                    if (z>=-5.0) then
!@                        C0 = 0.4+0.72*(5.0+z)
!@                    else
!@                        C0 = 0.4
!@                    end if
!@
!@                    mesh%mu_s (i,is,ief) = 0.18
!@                    mesh%mu_d (i,is,ief) = 0.12
!@                    mesh%Dc   (i,is,ief) = 0.30
!@                    mesh%C0   (i,is,ief) = C0
!@                    end if
!@!!!#endif
!@
!@!!!#ifdef TPV29
!@                    if ( trim(adjustl(problem)) .eq. 'tpv29' ) then
!@                    if(abs(z)<=17.0)then
!@                      omeg = 1.0
!@                    elseif(abs(z)<=22.0)then
!@                      omeg = (22.0 - abs(z))/5.0
!@                    else
!@                      omeg = 0.0
!@                    endif
!@                    Pf = 1.0*9.8*(-z) ! In MPa
!@                    szz = 2.670*9.8*(z)
!@                    syy = omeg*(1.025837*(szz+Pf)-Pf) + (1.0-omeg)*szz !\sigma_11
!@                    sxx = omeg*(0.974162*(szz+Pf)-Pf) + (1.0-omeg)*szz !\sigma_33
!@                    sxy =-omeg*(0.158649*(szz+Pf))                     !\sigma_13
!@                    sxy = -sxy
!@                    syz = 0
!@                    sxz = 0
!@
!@                    if (z>=-4.0) then
!@                        C0 = 0.4+0.2*(4.0+z)
!@                    else
!@                        C0 = 0.4
!@                    end if
!@
!@                    mesh%mu_s (i,is,ief) = 0.18
!@                    mesh%mu_d (i,is,ief) = 0.12
!@                    mesh%Dc   (i,is,ief) = 0.30
!@                    mesh%C0   (i,is,ief) = C0
!@                    end if
!@!!!#endif
!@
!@!!!#if defined (WENCHUAN)
!@                    if ( trim(adjustl(problem)) .eq. 'wenchuan' ) then
!@                    !if(abs(z)<=17.0)then
!@                    !  omeg = 1.0
!@                    !elseif(abs(z)<=22.0)then
!@                    !  omeg = (22.0 - abs(z))/5.0
!@                    !else
!@                    !  omeg = 0.0
!@                    !endif
!@                    !Pf = 1.0*9.8*(-z) ! In MPa
!@                    !szz = 2.670*9.8*(z)
!@                    !syy = omeg*(1.025837*(szz+Pf)-Pf) + (1.0-omeg)*szz !\sigma_11
!@                    !sxx = omeg*(0.974162*(szz+Pf)-Pf) + (1.0-omeg)*szz !\sigma_33
!@                    !sxy =-omeg*(0.158649*(szz+Pf))                     !\sigma_13
!@                    !sxy = -sxy
!@                    !syz = 0
!@                    !sxz = 0
!@                    depth = abs(z)
!@                    if(depth<10.e0)then
!@                      p_ratio=0.4*(1.0-cos(PI*depth/10.0e0)) +0.1 ! from 0.1 to 0.9
!@                    else
!@                      p_ratio=0.9
!@                    endif
!@                    if (depth<0.01) then
!@                      sigma_v = 0.01 / 8.0 * 12.0
!@                    elseif (depth<8.0) then
!@                      sigma_v = depth / 8.0 * 12.0
!@                    else
!@                      sigma_v = 12.0
!@                    endif
!@
!@                    T_local = 0.0
!@                    !T_local(3,3) = (1.0 - p_ratio)*2.670 * 9.8 * z
!@
!@                    ! sigma_H > sigma_h > sigma_v
!@                    T_local(3,3) = -sigma_v
!@                    T_local(2,2) = 1.1*T_local(3,3)
!@                    T_local(1,1) = 2.0*T_local(3,3)
!@
!@                    h1 = 4.0e0
!@                    h2 = 20.0e0
!@                    !if(y(i0,j,k)<305.0e3)then
!@                    !  h1 = 4.0e3;
!@                    !  h2 =15.0e3;
!@                    !else
!@                    !  h1 = 4.0e3;
!@                    !  h2 =25.0e3;
!@                    !endif
!@                    Dc0=0.3
!@                    if(depth<h1)then
!@                      Dc = (1.0+1.0*(1.0+cos(PI*(depth/h1))))*Dc0
!@                    elseif(depth>h2)then
!@                      Dc = ((depth-h2)/1.0e0*2.0+1.0)*Dc0
!@                    else
!@                      Dc = Dc0
!@                    endif
!@                    Dc = Dc0
!@
!@                    angle = 30.0*PI/180.0;
!@                    e_local(:,1) = (/ dcos(angle), dsin(angle), 0.0d0 /)
!@                    e_local(:,2) = (/-dsin(angle), dcos(angle), 0.0d0 /)
!@                    e_local(:,3) = (/ 0.0d0, 0.0d0,  1.0d0 /)
!@
!@                    e_global(:,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
!@                    e_global(:,2) = (/ 0.0d0, 1.0d0, 0.0d0 /)
!@                    e_global(:,3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
!@                    do m=1,3
!@                      do n=1,3
!@                        mat_L2G(m,n) = dot_product(e_local(:,m), e_global(:,n))
!@                      enddo
!@                    enddo
!@                    T_global = 0.0
!@                    do m=1,3
!@                      do n=1,3
!@                        T_global(1,1) = T_global(1,1) + mat_L2G(m,1)*mat_L2G(n,1)*T_local(m,n)
!@                        T_global(1,2) = T_global(1,2) + mat_L2G(m,1)*mat_L2G(n,2)*T_local(m,n)
!@                        T_global(1,3) = T_global(1,3) + mat_L2G(m,1)*mat_L2G(n,3)*T_local(m,n)
!@                        T_global(2,1) = T_global(2,1) + mat_L2G(m,2)*mat_L2G(n,1)*T_local(m,n)
!@                        T_global(2,2) = T_global(2,2) + mat_L2G(m,2)*mat_L2G(n,2)*T_local(m,n)
!@                        T_global(2,3) = T_global(2,3) + mat_L2G(m,2)*mat_L2G(n,3)*T_local(m,n)
!@                        T_global(3,1) = T_global(3,1) + mat_L2G(m,3)*mat_L2G(n,1)*T_local(m,n)
!@                        T_global(3,2) = T_global(3,2) + mat_L2G(m,3)*mat_L2G(n,2)*T_local(m,n)
!@                        T_global(3,3) = T_global(3,3) + mat_L2G(m,3)*mat_L2G(n,3)*T_local(m,n)
!@                      enddo
!@                    enddo
!@
!@                    sxx = T_global(1,1)
!@                    syy = T_global(2,2)
!@                    szz = T_global(3,3)
!@                    sxy = T_global(1,2)
!@                    sxz = T_global(1,3)
!@                    syz = T_global(2,3)
!@
!@                    if (z>=-4.0) then
!@                        C0 = 0.4+0.2*(4.0+z)
!@                    else
!@                        C0 = 0.4
!@                    end if
!@                    C0 = 0
!@
!@                    mesh%mu_s (i,is,ief) = 0.40
!@                    mesh%mu_d (i,is,ief) = 0.10
!@                    mesh%Dc   (i,is,ief) = Dc
!@                    mesh%C0   (i,is,ief) = C0
!@                    end if
!@!!!#endif
!@
!@                    Tx = sxx*vec_n(1)+sxy*vec_n(2)+sxz*vec_n(3)
!@                    Ty = sxy*vec_n(1)+syy*vec_n(2)+syz*vec_n(3)
!@                    Tz = sxz*vec_n(1)+syz*vec_n(2)+szz*vec_n(3)
!@                    call rotate_xyz2nml(vec_n,vec_m,vec_l,Tx,Ty,Tz,Tn,Tm,Tl)
!@
!@!!!#if defined (WENCHUAN)
!@                    if ( trim(adjustl(problem)) .eq. 'wenchuan' ) then
!@                    asp_size = 4e0
!@                    if ( abs(y-0e0)<=asp_size .and. abs(z+10e0)<=asp_size ) then
!@                        Tm = 1.0*abs(Tn) * (mesh%mu_s(i,is,ief)+0.01) * sign(1d0, Tm)
!@                        Tl = 1.0*abs(Tn) * (mesh%mu_s(i,is,ief)+0.01) * sign(1d0, Tl)
!@                    end if
!@                    end if
!@!!!#endif
!@
!@                    mesh%tau0n(i,is,ief) = Tn
!@                    mesh%tau0m(i,is,ief) = Tm
!@                    mesh%tau0l(i,is,ief) = Tl
!@
!@                    !if (abs(y+0e3)>9e3 .or. z<-9e3) then
!@                    !if (abs(y+0e0)>9e0 .or. abs(z+5)>4e0) then
!@                    !if (abs(yc+0e0)>=15e0 .or. zc<=-15e0) then
!@                    !    mesh%mu_s(i,is,ie) = 1e4
!@                    !    mesh%C0  (:,is,ie) = 1e9
!@                    !end if
!@
!@
!@                enddo
!@
!@            else
!@                mesh%mu_s(:,is,ief) = 1e4
!@                mesh%C0  (:,is,ief) = 1e9
!@            endif
!@
!@        enddo
!@    enddo
!@
!@    mesh%stress = mesh%tau0n
!@
!@    !if (.false.) then
!@    !write(filename,'(a,i6.6)') 'data/fault0mpi',myrank
!@    !open(100,file=trim(filename))
!@    !write(100,*) Nfp,mesh%nfault_elem
!@    !do ie = 1,mesh%nelem
!@    !do is = 1,Nfaces
!@    !if (mesh%bctype(is,ie) >= BC_FAULT) then
!@    !do i = 1,Nfp
!@    !    j = i+(is-1)*Nfp
!@    !    write(100,*) &
!@    !        mesh%vx(mesh%vmapM(i,is,ie)), &
!@    !        mesh%vy(mesh%vmapM(i,is,ie)), &
!@    !        mesh%vz(mesh%vmapM(i,is,ie)), &
!@    !        mesh%nx(j,ie), &
!@    !        mesh%ny(j,ie), &
!@    !        mesh%nz(j,ie), &
!@    !        mesh%tau0n(i,is,ie), &
!@    !        mesh%tau0m(i,is,ie), &
!@    !        mesh%tau0l(i,is,ie)
!@    !end do
!@    !end if
!@    !end do
!@    !end do
!@    !close(100)
!@    !end if
!@
!@
end subroutine fault_init

subroutine write_fault(mesh,it,myrank)
  implicit none
  type(meshVar),intent(in) :: mesh
  integer,intent(in) :: it,myrank
  character(len=128) :: filename
  integer :: i,is,ie
  write(filename,'(a,i6.6,a,i6.6)') 'data/fault_it',it,'_mpi',myrank
  open(100,file=trim(filename))
  write(100,*) Nfp,mesh%nfault_elem
  do ie = 1,mesh%Nelem
    do is = 1,Nfaces
      if (mesh%bctype(is,ie) >= BC_FAULT) then
        do i = 1,Nfp
          write(100,*) &
          mesh%vx(mesh%vmapM(i,is,ie)), &
          mesh%vy(mesh%vmapM(i,is,ie)), &
          mesh%vz(mesh%vmapM(i,is,ie)), &
          mesh%sliprate(i,is,ie), &
          mesh%stress  (i,is,ie), &
          mesh%ruptime (i,is,ie)
          !mesh%slip    (i,is,ie)
        end do
      end if
    end do
  end do
  close(100)
end subroutine

end module
