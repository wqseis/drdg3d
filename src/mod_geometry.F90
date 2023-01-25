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

module mod_geometry

  use mod_para,   only : RKIND,             &
                         Order, NGLL,       &
                         Np, Nfp, Nfaces,   &
                         EPS,               &
                         BC_FAULT
  use mod_nodes,  only : Nodes3D,           &
                         getGLL,            &
                         xyztorst
  use mod_jacobi, only : Vdm2D,             &
                         Vdm3D,             &
                         gradVdm3D,         &
                         invVdm2D,          &
                         invVdm3D,          &
                         Dmatrices3D
  use mod_types,  only : meshvar

  implicit none

contains

subroutine build_geometry(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: myrank
  real(kind=RKIND),dimension(Np) :: x,y,z,r,s,t
  real(kind=RKIND),dimension(Np) :: rx,sx,tx
  real(kind=RKIND),dimension(Np) :: ry,sy,ty
  real(kind=RKIND),dimension(Np) :: rz,sz,tz,jac
  real(kind=RKIND),dimension(3) :: n_n,n_m,n_l
  real(kind=RKIND) :: err,tmp1,tmp2,vxc,vyc,vzc
  integer :: Nelem
  integer :: ie,is,i,j,k,l,n
  integer :: va,vb,vc,vd
  integer :: iv(Np)

  Nelem = mesh%Nelem
  myrank = mesh%rank

  ! matrix
  allocate(mesh%Vdm(Np,Np))
  allocate(mesh%invVdm(Np,Np))
  allocate(mesh%Mass(Np,Np))
  allocate(mesh%Dr(Np,Np))
  allocate(mesh%Ds(Np,Np))
  allocate(mesh%Dt(Np,Np))
  allocate(mesh%Fmask(Nfp,Nfaces))
  allocate(mesh%LIFT(Np,Nfp*Nfaces))

  allocate(mesh%r(Np))
  allocate(mesh%s(Np))
  allocate(mesh%t(Np))
  allocate(mesh%vx(Np*Nelem))
  allocate(mesh%vy(Np*Nelem))
  allocate(mesh%vz(Np*Nelem))
  allocate(mesh%rx(Np*Nelem))
  allocate(mesh%ry(Np*Nelem))
  allocate(mesh%rz(Np*Nelem))
  allocate(mesh%sx(Np*Nelem))
  allocate(mesh%sy(Np*Nelem))
  allocate(mesh%sz(Np*Nelem))
  allocate(mesh%tx(Np*Nelem))
  allocate(mesh%ty(Np*Nelem))
  allocate(mesh%tz(Np*Nelem))
  allocate(mesh%jac(Np*Nelem))
  ! surface
  allocate(mesh%nx(Nfaces*Nfp,Nelem))
  allocate(mesh%ny(Nfaces*Nfp,Nelem))
  allocate(mesh%nz(Nfaces*Nfp,Nelem))
  allocate(mesh%mx(Nfaces*Nfp,Nelem))
  allocate(mesh%my(Nfaces*Nfp,Nelem))
  allocate(mesh%mz(Nfaces*Nfp,Nelem))
  allocate(mesh%lx(Nfaces*Nfp,Nelem))
  allocate(mesh%ly(Nfaces*Nfp,Nelem))
  allocate(mesh%lz(Nfaces*Nfp,Nelem))
  allocate(mesh%sJ(Nfaces*Nfp,Nelem))
  allocate(mesh%Fscale(Nfaces*Nfp,Nelem))
  ! media
  !allocate(mesh%vp (Nelem))
  !allocate(mesh%vs (Nelem))
  !allocate(mesh%rho(Nelem))
  !allocate(mesh%zp (Nelem))
  !allocate(mesh%zs (Nelem))
  !allocate(mesh%mu (Nelem))
  !allocate(mesh%lam(Nelem))
  ! connection
  allocate(mesh%vmapM(Nfp,Nfaces,Nelem))
  allocate(mesh%vmapP(Nfp,Nfaces,Nelem))

  call Nodes3D(x,y,z)
  call xyztorst(x,y,z,r,s,t)

  !print*,'r=',sngl(r)

  mesh%r = r
  mesh%s = s
  mesh%t = t

  !do j = 1,Np
  !  write(*,'(3f10.4)'),r(j),s(j),t(j)
  !end do
  
  do ie = 1,Nelem
    va = mesh%elem(1,ie)
    vb = mesh%elem(2,ie)
    vc = mesh%elem(3,ie)
    vd = mesh%elem(4,ie)
    do j = 1,Np
      mesh%vx(j+(ie-1)*Np) = 0.5d0*(             &
          -(1d0+r(j)+s(j)+t(j))*mesh%coord(1,va) &
          +(1d0+r(j)          )*mesh%coord(1,vb) &
          +(1d0+     s(j)     )*mesh%coord(1,vc) &
          +(1d0+          t(j))*mesh%coord(1,vd) )
      mesh%vy(j+(ie-1)*Np) = 0.5d0*(             &
          -(1d0+r(j)+s(j)+t(j))*mesh%coord(2,va) &
          +(1d0+r(j)          )*mesh%coord(2,vb) &
          +(1d0+     s(j)     )*mesh%coord(2,vc) &
          +(1d0+          t(j))*mesh%coord(2,vd) )
      mesh%vz(j+(ie-1)*Np) = 0.5d0*(             &
          -(1d0+r(j)+s(j)+t(j))*mesh%coord(3,va) &
          +(1d0+r(j)          )*mesh%coord(3,vb) &
          +(1d0+     s(j)     )*mesh%coord(3,vc) &
          +(1d0+          t(j))*mesh%coord(3,vd) )
    end do
  end do

  mesh%xmin = minval(mesh%vx)
  mesh%xmax = maxval(mesh%vx)
  mesh%ymin = minval(mesh%vy)
  mesh%ymax = maxval(mesh%vy)
  mesh%zmin = minval(mesh%vz)
  mesh%zmax = maxval(mesh%vz)

  !return
#ifdef DEBUG
  print*,'rank=',myrank,'Vdm3D'
#endif
  call Vdm3D(mesh%Vdm,r,s,t)
  call invVdm3D(mesh%Vdm,mesh%invVdm,0)
  mesh%Mass = matmul(transpose(mesh%invVdm),mesh%invVdm)
  call Dmatrices3D(mesh%Dr,mesh%Ds,mesh%Dt,r,s,t,mesh%Vdm)
  !do i = 1,Np
  !  print*,mesh%Dt(i,:)
  !end do

  !do ie = 1,1
  !  do j = 1,Np
  !    write(*,'(3f10.4)') &
  !      mesh%vx(j+(ie-1)*Np), &
  !      mesh%vy(j+(ie-1)*Np), &
  !      mesh%vz(j+(ie-1)*Np)
  !  end do
  !end do

#ifdef DEBUG
  print*,'rank=',myrank,'Fmask'
#endif
  ! find boundary nodes in the local element
  i=1;j=1;k=1;l=1;
  do n=1,Np
    if (abs(t(n)+1.0) < EPS) then
      mesh%Fmask(i,1)=n; i=i+1
    end if
    if (abs(s(n)+1.0) < EPS) then
      mesh%Fmask(j,2)=n; j=j+1
    end if
    if (abs(1.0+r(n)+s(n)+t(n)) < EPS) then
      mesh%Fmask(k,3)=n; k=k+1
    end if
    if (abs(r(n)+1.0) < EPS) then
      mesh%Fmask(l,4)=n; l=l+1
    end if
  end do

  !do i = 1,Nfp
  !  print*,mesh%Fmask(i,:)
  !end do

#ifdef DEBUG
  print*,'rank=',myrank,'Lift3D'
#endif
  call Lift3D(mesh%LIFT,mesh%Fmask,mesh%Vdm,r,s,t)

  !do i = 1,Np
  !  print*,mesh%LIFT(i,:)
  !end do
#ifdef DEBUG
  print*,'rank=',myrank,'GeometricFactors3D'
#endif
  do ie = 1,Nelem
    do j = 1,Np
      iv(j) = j+(ie-1)*Np
    end do
    call GeometricFactors3D(rx,sx,tx,ry,sy,ty,rz,sz,tz,jac,&
        mesh%vx(iv),mesh%vy(iv),mesh%vz(iv),&
        mesh%Dr,mesh%Ds,mesh%Dt)
    mesh%rx(iv) = rx
    mesh%sx(iv) = sx
    mesh%tx(iv) = tx
    mesh%ry(iv) = ry
    mesh%sy(iv) = sy
    mesh%ty(iv) = ty
    mesh%rz(iv) = rz
    mesh%sz(iv) = sz
    mesh%tz(iv) = tz
    mesh%jac(iv) = jac

    !if (ie==2) &
    !print*,sngl(ry);! stop 2
    !print*,abs(jac-jac(1));! stop 2
  end do

#ifdef DEBUG
  print*,'rank=',myrank,'Normals3D'
#endif
  ! get normal vectors
  do ie = 1,Nelem
    do j = 1,Np
      iv(j) = j+(ie-1)*Np
    end do
    call Normals3D(mesh%nx(:,ie),mesh%ny(:,ie),mesh%nz(:,ie),mesh%sJ(:,ie),&
        mesh%Dr,mesh%Ds,mesh%Dt,mesh%vx(iv),mesh%vy(iv),mesh%vz(iv),mesh%Fmask)

    if (ie==111 .and. .false.) then
    do j = 1,Nfp
      print*,&
        sngl(mesh%nz(j+0*Nfp,ie)),&
        sngl(mesh%nz(j+1*Nfp,ie)),&
        sngl(mesh%nz(j+2*Nfp,ie)),&
        sngl(mesh%nz(j+3*Nfp,ie))
    end do
    stop 
    end if
  end do

#ifdef DEBUG
  print*,'rank=',myrank,'tangential vectors'
#endif
  ! tangential vectors
  do ie = 1,Nelem
    do is = 1,Nfaces
      do i = 1,Nfp
        n_n = (/ &
            mesh%nx(i+(is-1)*Nfp,ie), &
            mesh%ny(i+(is-1)*Nfp,ie), &
            mesh%nz(i+(is-1)*Nfp,ie) /)
        call gen_orth_vectors(n_n, n_m, n_l)
        mesh%mx(i+(is-1)*Nfp,ie) = n_m(1)
        mesh%my(i+(is-1)*Nfp,ie) = n_m(2)
        mesh%mz(i+(is-1)*Nfp,ie) = n_m(3)
        mesh%lx(i+(is-1)*Nfp,ie) = n_l(1)
        mesh%ly(i+(is-1)*Nfp,ie) = n_l(2)
        mesh%lz(i+(is-1)*Nfp,ie) = n_l(3)
      end do
    end do
  end do

#ifdef DEBUG
  print*,'rank=',myrank,'Fscale vectors'
#endif
  ! make scaling vector
  do ie = 1,Nelem
    do j = 1,Np
      iv(j) = j+(ie-1)*Np
    end do
    do is=1,Nfaces
      do j=1,Nfp
        k=iv(mesh%Fmask(j,is))
        mesh%Fscale(j+(is-1)*Nfp,ie)=mesh%sJ(j+(is-1)*Nfp,ie)/mesh%Jac(k)
      end do
    end do

    if (ie==123 .and. .false.) then
      do j = 1,Nfp
        print*,&
          sngl(mesh%Fscale(j+0*Nfp,ie)),&
          sngl(mesh%Fscale(j+1*Nfp,ie)),&
          sngl(mesh%Fscale(j+2*Nfp,ie)),&
          sngl(mesh%Fscale(j+3*Nfp,ie))
      end do
      stop
    end if
  end do

#ifdef DEBUG
  print*,'rank=',myrank,'build_maps3D'
#endif
  call build_maps3D(mesh)

#ifdef DEBUG
  print*,'rank=',myrank,'checking neighbours'
#endif
  ! checking mappings
  do ie = 1,Nelem
    do is = 1,Nfaces
      if(mesh%neigh(is,ie) >= 0) then ! do not check mpi
        err =  &
            sum(abs(mesh%vx(mesh%vmapM(:,is,ie))-mesh%vx(mesh%vmapP(:,is,ie)))) + &
            sum(abs(mesh%vy(mesh%vmapM(:,is,ie))-mesh%vy(mesh%vmapP(:,is,ie)))) + &
            sum(abs(mesh%vz(mesh%vmapM(:,is,ie))-mesh%vz(mesh%vmapP(:,is,ie)))) 
        err = err/Nfp
        !if(err>1e-6) then
        if(err>1e-3) then
          print*,'check mapping error ',ie,mesh%neigh(is,ie),is,sngl(err)
          stop 17
        end if
      end if
    end do
  end do

  !call write_mesh(mesh)
  !call write_normals(mesh,0)

!!!!  mesh%vp (:) = 6.000
!!!!  mesh%vs (:) = 3.464
!!!!  mesh%rho(:) = 2.670
!!!!#ifdef TPV6
!!!!  do ie = 1,Nelem
!!!!    do j = 1,Np
!!!!      iv(j) = j+(ie-1)*Np
!!!!    end do
!!!!    vxc = sum(mesh%vx(iv))/dble(Np)
!!!!    if(vxc<0) then
!!!!      mesh%vp (ie) = 6.0/1.6
!!!!      mesh%vs (ie) = 3.464/1.6
!!!!      mesh%rho(ie) = 2.670/1.2
!!!!    end if
!!!!  end do
!!!!#endif

  mesh%zp = mesh%rho*mesh%vp
  mesh%zs = mesh%rho*mesh%vs
  mesh%mu = mesh%rho*mesh%vs**2
  mesh%lam = mesh%rho*mesh%vp**2-2.0*mesh%mu

#ifdef DEBUG
  print*,'rank=',myrank,'setDT'
#endif
  call setDT(mesh)

#ifdef DEBUG
  print*,'rank=',myrank,'set bc'
#endif
  ! set bctype
  do ie = 1,Nelem
    do is = 1,Nfaces
      !err =sum(abs(mesh%vz(mesh%vmapM(:,is,ie))-0.5))/Nfp 
      err =sum(abs(mesh%vx(mesh%vmapM(:,is,ie))-0.0))/Nfp 
      err =sum(abs(mesh%vy(mesh%vmapM(:,is,ie))-0.0))/Nfp 
      vxc =sum((mesh%vx(mesh%vmapM(:,is,ie))-0.0))/Nfp 
      vyc =sum((mesh%vy(mesh%vmapM(:,is,ie))-0.0))/Nfp 
      vzc =sum((mesh%vz(mesh%vmapM(:,is,ie))-0.0))/Nfp 
!!!#if defined (TPV5) || defined (TPV6)
!!!      !if(err<1e-3) then
!!!      if(abs(vxc)<1e-3 .and. &
!!!          vyc>=-15 .and. vyc<=15 .and. &
!!!          vzc>=-15) then
!!!        mesh%bctype(is,ie) = BC_FAULT
!!!      end if
!!!#endif

!!!#ifdef TPV24
!!!      ! main fault
!!!      if(abs(vxc)<1e-3 .and. &
!!!          vyc>=-16 .and. vyc<=12 .and. &
!!!          vzc>=-15) then
!!!        mesh%bctype(is,ie) = BC_FAULT
!!!      end if
!!!      ! branch fault
!!!      if( abs(vyc-vxc*dsqrt(3d0)) < 1e-3 .and. & 
!!!          vxc >= 0 .and. vxc <= 6 .and. &
!!!          vzc >= -15) then
!!!        mesh%bctype(is,ie) = BC_FAULT
!!!      end if
!!!#endif
!!!
!!!#ifdef TPV26
!!!      if(abs(vxc)<1e-3 .and. &
!!!          vyc>=-20 .and. vyc<=20 .and. &
!!!          vzc>=-20) then
!!!        mesh%bctype(is,ie) = BC_FAULT
!!!      end if
!!!#endif

      ! free surface
      !if(abs(vzc)<1e-3) then
      !  mesh%bctype(is,ie) = BC_FREE
      !end if
    end do
  end do

  !stop
  if(.false.) then
  ! check normal vector
  do ie = 1,Nelem
    do is = 1,Nfaces
      if(mesh%bctype(is,ie) >= BC_FAULT) then
        !print*,'ie=',ie,'nx=',mesh%nx(1+(is-1)*Nfp:is*Nfp,ie)
        !print*,'ie=',ie,'nx=',mesh%Fscale(1+(is-1)*Nfp:is*Nfp,ie)
        !print*,'ie=',ie,'rx=',mesh%rx(mesh%vmapM(:,is,ie))
        tmp1=mesh%jac(mesh%vmapM(1,is,ie))
        tmp2=mesh%jac(mesh%vmapP(1,is,ie))
        if(tmp1<tmp2) then
          tmp1=tmp2/tmp1
        else
          tmp1=tmp1/tmp2
        endif
        !print*,'ie=',ie,'neig=',mesh%neigh(is,ie),'jac=',tmp1
        do i = 1,Nfp
          if(mesh%nx(i+(is-1)*Nfp,ie) .gt. 0) then
            tmp1=mesh%jac(mesh%vmapM(i,is,ie))
            tmp2=mesh%jac(mesh%vmapP(i,is,ie))
            if(tmp1<tmp2) then
              tmp1=tmp2/tmp1
            else
              tmp1=tmp1/tmp2
            endif
            print*,&
            mesh%vx(mesh%vmapM(i,is,ie)),&
            mesh%vy(mesh%vmapM(i,is,ie)),&
            mesh%vz(mesh%vmapM(i,is,ie)),&
            tmp1
          end if
        end do
      end if
    end do
  end do
  stop 2
  end if
  !stop

end subroutine

subroutine build_maps3D(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: N,i,j,k,Nelem,is,ie,dire,face,neigh
  !integer :: FtoV(4,3)!,Fmask(Nfp,Nfaces)
  integer, allocatable :: tri(:,:),tri1(:),tri2(:),tri3(:)
  integer, allocatable :: ttri1(:),ttri2(:),ttri3(:)
  integer, allocatable :: tris(:,:)
  !data((FtoV(i,j),i=1,4),j=1,3)/1,1,2,1,2,2,3,3,3,4,4,4/
  !FtoV = reshape((/1,1,2,1,2,2,3,3,3,4,4,4/),shape(FtoV))

  Nelem = mesh%Nelem
  N=Order+1
  allocate(tri(N,N))
  k=1
  do i=1,N
    do j=1,N+1-i
      tri(j,i)=k; k=k+1
    enddo
  enddo
  k=(N*(N+1))/2
  allocate( tri1(k));allocate( tri2(k));allocate( tri3(k))
  allocate(ttri1(k));allocate(ttri2(k));allocate(ttri3(k))
  tri1=0;tri2=0;tri3=0;ttri1=0;ttri2=0;ttri3=0;k=1
  do i=1,N
      do j=1,N+1-i
          tri1(tri(i,j))=k
          tri3(tri(j,N+2-i-j))=k
          tri2(tri(N+2-i-j,i))=k
          ttri1(tri(j,i))=k
          ttri3(tri(N+2-i-j,j))=k
          ttri2(tri(i,N+2-i-j))=k
          k=k+1
      enddo
  enddo

  k = Nfp
  allocate(tris(k,6))
  tris(:,1) = tri1
  tris(:,2) = tri2
  tris(:,3) = tri3
  tris(:,4) = ttri1
  tris(:,5) = ttri2
  tris(:,6) = ttri3

  allocate(mesh%flipped_index(Nfp,6))
  mesh%flipped_index = tris

  do ie = 1,Nelem
    do is = 1,Nfaces
      neigh = mesh%neigh(is,ie)
      face = mesh%face(is,ie)
      dire = mesh%direction(is,ie)
      mesh%vmapM(:,is,ie) = mesh%Fmask(:,is) + (ie-1)*Np
      if(neigh==0 .or. neigh==ie) then
        mesh%vmapP(:,is,ie) = mesh%vmapM(:,is,ie)
      else if (neigh>0) then
        mesh%vmapP(:,is,ie) = mesh%Fmask(tris(:,dire),face) + (neigh-1)*Np
      end if
    end do
  end do

  deallocate(tri);
  deallocate(tri1);
  deallocate(tri2);
  deallocate(tri3);
  deallocate(ttri1);
  deallocate(ttri2);
  deallocate(ttri3);
  deallocate(tris);

end subroutine

subroutine Lift3D(LIFT,Fmask,Vdm,r,s,t)
  implicit none
  real(kind=RKIND), dimension(:), intent(in) :: r,s,t
  !real(kind=RKIND), dimension(Nfp,Nfp) :: v,m,inv
  real(kind=RKIND), dimension(:,:), allocatable :: v,m,inv
  real(kind=RKIND), dimension(:,:), intent(in) :: vdm
  integer, dimension(:,:), intent(in) :: Fmask
  !real(kind=RKIND), dimension(Np,Nfaces*Nfp), intent(out) :: LIFT
  real(kind=RKIND), dimension(:,:), intent(out) :: LIFT
  integer :: iface,i
  !real(kind=RKIND), dimension(Np,Nfaces*Nfp) :: emat,temp
  real(kind=RKIND), dimension(:,:),allocatable :: emat,temp
  real(kind=RKIND), dimension(Nfp) :: faceR, faceS
  integer, dimension(Nfp) :: idr,idc,id

  allocate(v(Nfp,Nfp))
  allocate(m(Nfp,Nfp))
  allocate(inv(Nfp,Nfp))
  allocate(emat(Np,Nfp*Nfaces))
  allocate(temp(Np,Nfp*Nfaces))

  emat(:,:) = 0

  do iface=1,Nfaces
    v=0; inv=0; m=0
    if (iface == 1) then
      id=Fmask(:,1); faceR=r(id); faceS=s(id)
    else if (iface == 2) then
      id=Fmask(:,2); faceR=r(id); faceS=t(id)
    else if (iface == 3) then
      id=Fmask(:,3); faceR=s(id); faceS=t(id)
    else if (iface == 4) then
      id=Fmask(:,4); faceR=s(id); faceS=t(id)
    end if

    call Vdm2D(v,faceR,faceS)
    call invVdm2D(v,inv,0,Nfp)

    m = matmul(transpose(inv),inv)

    idr=Fmask(:,iface)
    do i=1,Nfp
      idc(i)=((iface-1)*Nfp+1)+(i-1)
    end do
    emat(idr,idc) = emat(idr,idc)+m
  end do

  temp = matmul(transpose(vdm),emat)
  LIFT = matmul(Vdm,temp)

  deallocate(v,inv,m,emat,temp)
end subroutine Lift3D

subroutine GeometricFactors3D(rx,sx,tx,ry,sy,ty,rz,sz,tz,J,x,y,z,Dr,Ds,Dt)
  ! INPUT: differentation matrices Dr, Ds, Dt, global coordinates of points x,y,z
  ! DOES: compute the metric elements for local mappings of an element
  ! OUTPUT: geometric factors rx,sx,tx,ry,sy,ty,rz,sz,tz; Jacobian J
  implicit none
  real(kind=RKIND), dimension(:), intent(in) :: x,y,z
  real(kind=RKIND), dimension(:,:), intent(in) :: Dr,Ds,Dt
  real(kind=RKIND), dimension(size(x)) :: rx,sx,tx,ry,sy,ty,rz,sz,tz,J
  real(kind=RKIND), dimension(size(x)) :: xr,xs,xt,yr,ys,yt,zr,zs,zt
  !integer :: i

  xr = matmul(Dr,x)               ! derivative of with respect to r
  xs = matmul(Ds,x)
  xt = matmul(Dt,x)
  yr = matmul(Dr,y)
  ys = matmul(Ds,y)
  yt = matmul(Dt,y)
  zr = matmul(Dr,z)
  zs = matmul(Ds,z)
  zt = matmul(Dt,z)

  J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt)  ! get Jacobian

  !do i=1,size(J)
  !   if ( abs(J(i) - 0.0) < EPS) then
  !      write(*,*)  i, J(i)
  !      stop "Jacobian is zero in geometric factors"
  !   end if
  !end do

  rx =  (ys*zt - zs*yt)/J     ! inverse mapping (derivative of r with respect to x)
  ry = -(xs*zt - zs*xt)/J
  rz =  (xs*yt - ys*xt)/J

  sx = -(yr*zt - zr*yt)/J
  sy =  (xr*zt - zr*xt)/J
  sz = -(xr*yt - yr*xt)/J

  tx =  (yr*zs - zr*ys)/J
  ty = -(xr*zs - zr*xs)/J
  tz =  (xr*ys - yr*xs)/J
end subroutine GeometricFactors3D

subroutine Normals3D(nx,ny,nz,sJ,Dr,Ds,Dt,x,y,z,Fmask)
  ! INPUT: differentation matrices Dr,Ds,Dt, global coordinates x,y,z, surface mask fmask
  ! DOES: computes normals and their absolute value
  ! OUTPUT: normals nx,ny,nz and absolute value sj
  implicit none
  real(kind=RKIND), dimension(:,:), intent(in) :: Dr,Ds,Dt
  real(kind=RKIND), dimension(:), intent(in) :: x,y,z
  integer, dimension(:,:), intent(in) :: Fmask
  real(kind=RKIND), dimension(Nfaces*Nfp), intent(out) :: nx,ny,nz,sJ
  real(kind=RKIND), dimension(Np) :: xr,yr,zr,xs,ys,zs,xt,yt,zt,J
  real(kind=RKIND), dimension(Np) :: rx,sx,tx,ry,sy,ty,rz,sz,tz
  real(kind=RKIND), dimension(Nfaces*Nfp) :: frx,fry,frz,fsx,fsy,fsz,ftx,fty,ftz

  integer, dimension(Nfaces*Nfp) :: Fmaskv
  integer, dimension(Nfp) :: fid1,fid2,fid3,fid4
  integer :: i,k,l

  k=1
  do l=1,4       ! for all surfaces of tet
    do i=1,Nfp  ! for all points on surface
      Fmaskv(k)=Fmask(i,l) ! make vector with points on surface and 0 if specific point is not on surface
      if (l==1) then
        fid1(i)=k     ! vector with 4 blocks, first number 1 to Nfp, rest zeros
      else if (l==2) then
        fid2(i)=k     ! vector with 4 blocks, first zeros, than numbers Nfp+1 to 2*Nfp, than zeros again
      else if (l==3) then
        fid3(i)=k     ! like above (in third block)
      else if (l==4) then
        fid4(i)=k     ! like above (in fourth block)
      end if
      k=k+1
    end do
  end do

  xr = matmul(Dr,x)
  xs = matmul(Ds,x)
  xt = matmul(Dt,x)
  yr = matmul(Dr,y)
  ys = matmul(Ds,y)
  yt = matmul(Dt,y)
  zr = matmul(Dr,z)
  zs = matmul(Ds,z)
  zt = matmul(Dt,z)

  J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt)

  !do i=1,size(J)
  !   !if ( abs(J(i) - 0.0) < EPS) then
  !   if ( abs(J(i) - 0.0) < 1e-16) then
  !      write(*,*)  i, J(i)
  !      stop "jacobian zero in normals"
  !   end if
  !end do

  rx =  (ys*zt - zs*yt)/J
  ry = -(xs*zt - zs*xt)/J
  rz =  (xs*yt - ys*xt)/J

  sx = -(yr*zt - zr*yt)/J
  sy =  (xr*zt - zr*xt)/J
  sz = -(xr*yt - yr*xt)/J

  tx =  (yr*zs - zr*ys)/J
  ty = -(xr*zs - zr*xs)/J
  tz =  (xr*ys - yr*xs)/J

  ! interpolate geometric factors to face nodes
  frx = rx(Fmaskv)
  fsx = sx(Fmaskv)
  ftx = tx(Fmaskv)
  fry = ry(Fmaskv)
  fsy = sy(Fmaskv)
  fty = ty(Fmaskv)
  frz = rz(Fmaskv)
  fsz = sz(Fmaskv)
  ftz = tz(Fmaskv)

  ! build normals (use derivatives at surface points as normals)

  ! face 1
  nx(fid1) = -ftx(fid1)
  ny(fid1) = -fty(fid1)
  nz(fid1) = -ftz(fid1)
  ! face 2
  nx(fid2) = -fsx(fid2)
  ny(fid2) = -fsy(fid2)
  nz(fid2) = -fsz(fid2)
  ! face 3
  nx(fid3) = frx(fid3)+fsx(fid3)+ftx(fid3)
  ny(fid3) = fry(fid3)+fsy(fid3)+fty(fid3)
  nz(fid3) = frz(fid3)+fsz(fid3)+ftz(fid3)
  !face 4
  nx(fid4) = -frx(fid4)
  ny(fid4) = -fry(fid4)
  nz(fid4) = -frz(fid4)

  ! normalize
  !sj = dsqrt(nx*nx + ny*ny + nz*nz)
  sj = sqrt(nx*nx + ny*ny + nz*nz)
  nx = nx/sJ
  ny = ny/sJ
  nz = nz/sJ

  sJ = sJ * J(Fmaskv(:))
end subroutine Normals3D

subroutine gen_orth_vectors(n, m, l)
  implicit none

  real(kind=RKIND),dimension(:),intent(in) :: n
  real(kind=RKIND),dimension(:),intent(inout) :: l, m
  real(kind=RKIND) :: tol, diff_norm1, diff_norm2

  tol = 1.0e-12

  m(1) = 0d0
  m(2) = 1d0
  m(3) = 0d0

  !if (abs(sqrt(n(1)**2+n(2)**2+n(3)**)-1d0) .ge. tol) stop 'input vector must be a unit vector'

  !diff_norm1 = dsqrt((n(1)-m(1))**2 + (n(2)-m(2))**2 + (n(3)-m(3))**2)
  !diff_norm2 = dsqrt((n(1)+m(1))**2 + (n(2)+m(2))**2 + (n(3)+m(3))**2)
  diff_norm1 = sqrt((n(1)-m(1))**2 + (n(2)-m(2))**2 + (n(3)-m(3))**2)
  diff_norm2 = sqrt((n(1)+m(1))**2 + (n(2)+m(2))**2 + (n(3)+m(3))**2)

  if (diff_norm1 .ge. tol .and. diff_norm2 .ge. tol) then
    call Gram_Schmidt(n, m)
  else
    m(1) = 0d0
    m(2) = 0d0
    m(3) = 1d0
    !diff_norm1 = dsqrt((n(1)-m(1))**2 + (n(2)-m(2))**2 + (n(3)-m(3))**2)
    !diff_norm2 = dsqrt((n(1)+m(1))**2 + (n(2)+m(2))**2 + (n(3)+m(3))**2)
    diff_norm1 = sqrt((n(1)-m(1))**2 + (n(2)-m(2))**2 + (n(3)-m(3))**2)
    diff_norm2 = sqrt((n(1)+m(1))**2 + (n(2)+m(2))**2 + (n(3)+m(3))**2)
    if (diff_norm1 .ge. tol .and. diff_norm2 .ge. tol) then
      call Gram_Schmidt(n, m)
    else
      m(1) = 1d0
      m(2) = 0d0
      m(3) = 0d0
      call Gram_Schmidt(n, m )
    end if
  end if

  l(1) = n(2)*m(3)-n(3)*m(2)
  l(2) = n(3)*m(1)-n(1)*m(3)
  l(3) = n(1)*m(2)-n(2)*m(1)
end subroutine gen_orth_vectors

subroutine setDT(mesh)
  implicit none
  type(meshvar) :: mesh

  real(kind=RKIND) :: hmin,hmax,h,minGLL,xgll(NGLL)
  integer :: ie,i
  integer :: iv(Np)

  hmin = 1e30
  hmax = 0

  ! find min length

  do ie = 1,mesh%Nelem
    do i = 1,Np
      iv(i) = i+(ie-1)*Np
    end do
    h = minval(mesh%Jac(iv))/minval(mesh%sJ(:,ie))
    hmin = min(h,hmin)
    h = maxval(mesh%Jac(iv))/maxval(mesh%sJ(:,ie))
    hmax = max(h,hmax)
  end do

  !print*,'h=',hmin,'~',hmax
  mesh%rmin = hmin
  mesh%rmax = hmax

  xgll = getGLL()
  minGLL = xgll(1)-xgll(2)
  !print*,'xgll=',sngl(xgll)
  !print*,'minGLL=',minGLL

  mesh%dtfactor = 2.0/3.0*minGLL*hmin/maxval(mesh%vp)

  !print*,'dtfactor=',mesh%dtfactor

end subroutine

subroutine Gram_Schmidt(y, z)
  implicit none
  ! Gram Schmidt orthonormalization
  real(kind=rkind), dimension(:), intent(in) :: y
  real(kind=rkind), dimension(:), intent(inout) :: z
  real(kind=rkind) :: a_yz

  a_yz = y(1)*z(1) + y(2)*z(2) + y(3)*z(3)
  z = z - a_yz*y
  !z = z/dsqrt(z(1)**2 + z(2)**2 + z(3)**2)
  z = z/sqrt(z(1)**2 + z(2)**2 + z(3)**2)
end subroutine Gram_Schmidt

subroutine write_mesh(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: i,Nelem
  Nelem = mesh%Nelem
  open(100,file='vertices.txt')
  !write(100,*) Nelem*Np
  write(100,*) Np,Nelem
  do i = 1,Np*Nelem
    write(100,*) mesh%vx(i),mesh%vy(i),mesh%vz(i)
  enddo
  close(100)
end subroutine

subroutine write_normals(mesh,n)
  implicit none
  type(meshvar),intent(in) :: mesh
  integer,intent(in) :: n
  integer :: i,ie,is
  character(len=80) :: filename
  write(filename,'(a,i6.6)') 'data/normals',n
  !write(filename,'(a)') 'normals.txt'
  !open(100,file='data/vertices.txt')
  open(100,file=trim(filename))
  write(100,*) Nfp,Nfaces,mesh%nelem
  do ie = 1,mesh%nelem
    do is = 1,Nfaces
      do i = 1,Nfp
        write(100,*) &
        mesh%vx(mesh%vmapM(i,is,ie)),&
        mesh%vy(mesh%vmapM(i,is,ie)),&
        mesh%vz(mesh%vmapM(i,is,ie)),&
        mesh%nx(i+(is-1)*Nfp,ie),&
        mesh%ny(i+(is-1)*Nfp,ie),&
        mesh%nz(i+(is-1)*Nfp,ie),&
        mesh%mx(i+(is-1)*Nfp,ie),&
        mesh%my(i+(is-1)*Nfp,ie),&
        mesh%mz(i+(is-1)*Nfp,ie),&
        mesh%lx(i+(is-1)*Nfp,ie),&
        mesh%ly(i+(is-1)*Nfp,ie),&
        mesh%lz(i+(is-1)*Nfp,ie)
      end do
    end do
  end do
  close(100)
end subroutine

end module
