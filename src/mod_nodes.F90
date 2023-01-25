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

module mod_nodes

  use mod_para, only : RKIND,           &
                       CUSTOM_REAL,     &
                       Order, NGLL, Np, &
                       EPS,             &
                       PI
  !use mod_gll
  !use mod_warpfactor
  use mod_matrix

  implicit none

contains

subroutine Nodes3D(x,y,z)
  implicit none

  real(kind=RKIND),dimension(15),parameter :: alphastore = (/ &
          0d0,0d0,0d0,0.1002d0,1.1332d0,1.5608d0, &
          1.3413d0,1.2577d0,1.1603d0,1.10153d0,0.6080d0, &
          0.4523d0,0.8856d0,0.8717d0,0.9655d0/)
  !real(kind=RKIND), intent(in) :: alpha
  real(kind=RKIND) :: alpha
  real(kind=RKIND), dimension(:), intent(out) :: x,y,z
  real(kind=RKIND), dimension(Np) :: L1,L2,L3,L4,r,s,t
  real(kind=RKIND), dimension(Np) :: La,Lb,Lc,Ld
  real(kind=RKIND), dimension(Np) :: warp1,warp2
  real(kind=RKIND), dimension(Np) :: blend,denom
  integer, dimension(3) :: bool
  integer :: p=order
  real(kind=RKIND) :: tol=1e-10
  integer :: sk,n,m,q,i,j
  real(kind=RKIND), dimension(3) :: v1,v2,v3,v4
  real(kind=RKIND), dimension(Np,3) :: XYZ,shift
  real(kind=RKIND), dimension(4,3) :: t1,t2

  if(Order .le. 15)then
    alpha=alphastore(Order)
  else
    alpha=1.0d0
  endif

  !make eqidistant point in tet
  sk = 1
  do n=1,p+1
    do m=1,p+2-n
      do q=1,p+3-n-m
        r(sk) = -1d0 + (q-1d0)*2d0/p
        s(sk) = -1d0 + (m-1d0)*2d0/p
        t(sk) = -1d0 + (n-1d0)*2d0/p
        sk = sk+1
      end do
    end do
  end do

  L1(:) = (1d0+t(:))/2d0
  L2(:) = (1d0+s(:))/2d0
  L3(:) = -(1d0+r(:)+s(:)+t(:))/2d0
  L4(:) = (1d0+r(:))/2d0

  ! set vertices of the tetra
  v1 = (/-1d0, -1d0/dsqrt(3d0), -1d0/dsqrt(6d0) /)
  v2 = (/ 1d0, -1d0/dsqrt(3d0), -1d0/dsqrt(6d0) /)
  v3 = (/ 0d0,  2d0/dsqrt(3d0), -1d0/dsqrt(6d0) /)
  v4 = (/ 0d0,             0d0,  3d0/dsqrt(6d0) /)

  ! set orthogonal axis tengents on the faces
  t1(1,:) = v2(:) - v1(:)
  t1(2,:) = v2(:) - v1(:)
  t1(3,:) = v3(:) - v2(:)
  t1(4,:) = v3(:) - v1(:)
  t2(1,:) = v3(:) - 0.5d0*(v1(:)+v2(:))
  t2(2,:) = v4(:) - 0.5d0*(v1(:)+v2(:))
  t2(3,:) = v4(:) - 0.5d0*(v2(:)+v3(:))
  t2(4,:) = v4(:) - 0.5d0*(v1(:)+v3(:))

  ! normalise tangents
  do i=1,4
    t1(i,:) = t1(i,:)/sqrt(t1(i,1)**2 + t1(i,2)**2 + t1(i,3)**2)
    t2(i,:) = t2(i,:)/sqrt(t2(i,1)**2 + t2(i,2)**2 + t2(i,3)**2)
  end do

  ! warp and blend for each face
  shift(:,:) = 0d0

  do i=1,Np
    do j=1,3
      XYZ(i,j)=L3(i)*v1(j) + L4(i)*v2(j) + L2(i)*v3(j) + L1(i)*v4(j)
    end do
  end do

  do i=1,4
    if     (i == 1) then; La = L1; Lb = L2; Lc = L3; Ld = L4
    else if(i == 2) then; La = L2; Lb = L1; Lc = L3; Ld = L4
    else if(i == 3) then; La = L3; Lb = L1; Lc = L4; Ld = L2
    else                ; La = L4; Lb = L1; Lc = L3; Ld = L2
    end if
    ! compute warp tangential to face
    call shiftwarp(alpha,warp1,warp2,Lb,Lc,Ld)

    ! volume blending
    blend(:) = Lb(:)*Lc(:)*Ld(:)
    denom(:) = (Lb(:)+0.5d0*La(:))*(Lc(:)+0.5d0*La(:))*(Ld(:)+0.5d0*La(:))
    do j=1,Np
      if (denom(j) > tol) then
        blend(j) = ( 1d0 + (alpha*La(j))**2 ) * blend(j)/denom(j)
      end if
    end do

    !compute warp and shift
    do j=1,3
      shift(:,j) = shift(:,j) + (blend(:) *warp1) * t1(i,j) + (blend(:) *warp2) * t2(i,j)
    end do

    ! fix face warp
    do j=1,Np
      bool(:)=0
      if (Lb(j) > tol) bool(1)=1
      if (Lc(j) > tol) bool(2)=1
      if (Ld(j) > tol) bool(2)=1
      if ( (La(j) < tol) .AND. ( sum(bool)<3 )) then
        do n=1,3
          shift(j,n) = warp1(j) *t1(i,n) + warp2(j) *t2(i,n)
        end do
      end if
    end do
  end do
  XYZ = XYZ + shift
  x=XYZ(:,1)
  y=XYZ(:,2)
  z=XYZ(:,3)

end subroutine

function warpfactor(xnodes,xout)
  integer :: p
  !real(kind=RKIND), dimension(NGLL) :: xnodes,xeq
  !real(kind=RKIND), dimension(Np) :: warpfactor,d, warp ,xout
  real(kind=RKIND),allocatable,dimension(:) :: xnodes,xeq
  real(kind=RKIND),allocatable,dimension(:) :: warpfactor,d, warp ,xout
  real(kind=RKIND) :: dx,prod,prod1,prod2,e,sf
  logical :: zerof
  integer :: i,j,k, zeroff

  !allocate(xnodes(NGLL))
  allocate(xeq(NGLL))
  allocate(warpfactor(Np))
  allocate(d(Np))
  allocate(warp(Np))
  !allocate(xout(Np))

  p=NGLL-1

! create equispaced vector in same ordering as xnodes
  dx=2./(NGLL-1)
  j=NGLL
  do i=1,NGLL
    xeq(i)=((j-1.0)*dx)-1.0
    j=j-1
  end do
  d(:)=0.
  do k=1,Np
    do i=1,NGLL
      e=(xnodes(i)-xeq(i))
      prod1=1.
      prod2=1.
      do j=1,NGLL
        if (i/=j) then
          prod1=prod1*(xout(k)-xeq(j))
          prod2=prod2*(xeq(i)-xeq(j))
        end if
      end do
      prod=e*(prod1/prod2)
      d(k)=d(k)+prod
    end do
    warp(k)=d(k)

! scale warp factor
    zerof = abs(xout(k)) < (1.0-EPS)
    if (zerof) then
      zeroff=1.0
    else
      zeroff=0.0
    end if
    sf = 1.0-(zeroff*xout(k))**2
    warpfactor(k)=warp(k)/sf + warp(k)*(zeroff-1.0)
  end do
  deallocate(xeq,d,warp)
end function warpfactor

subroutine shiftwarp(alpha,x,y,L1,L2,L3)
  implicit none
  real(kind=RKIND), intent(in) :: alpha
  real(kind=RKIND), dimension(:), intent(out) :: x,y
  real(kind=RKIND), dimension(Np) :: L1,L2,L3, blend1, blend2, blend3
  real(kind=RKIND), dimension(Np) :: warpfactor1,warpfactor2,warpfactor3
  real(kind=RKIND), dimension(Np) :: warp1,warp2,warp3
  !real(kind=RKIND), dimension(NGLL) :: xgll
  real(kind=RKIND),allocatable,dimension(:) :: xgll
  integer :: i,p
  real(kind=RKIND),allocatable,dimension(:) :: dL

  p=NGLL-1

  allocate(xgll(NGLL))
  allocate(dL(Np))

  ! get GLL nodes
  xgll=getGLL()

  ! compute blending
  blend1(:) = 4.0*L2(:)*L3(:)
  blend2(:) = 4.0*L1(:)*L3(:)
  blend3(:) = 4.0*L1(:)*L2(:)

  !get warpfactors
  do i = 1,Np
    dL(i) = L3(i)-L2(i)
  end do
  warpfactor1 = warpfactor(xgll,dL)
  do i = 1,Np
    dL(i) = L1(i)-L3(i)
  end do
  warpfactor2 = warpfactor(xgll,dL)
  do i = 1,Np
    dL(i) = L2(i)-L1(i)
  end do
  warpfactor3 = warpfactor(xgll,dL)
  !warpfactor1 = warpfactor(xgll, L3(:)-L2(:))
  !warpfactor2 = warpfactor(xgll, L1(:)-L3(:))
  !warpfactor3 = warpfactor(xgll, L2(:)-L1(:))

  ! combine blend and warp
  warp1(:) = blend1(:)*warpfactor1(:)*(1.0+(alpha*L1(:))**2)
  warp2(:) = blend2(:)*warpfactor2(:)*(1.0+(alpha*L2(:))**2)
  warp3(:) = blend3(:)*warpfactor3(:)*(1.0+(alpha*L3(:))**2)

  ! accumulate deformations associated with each edge
  x(:) =  1.0*warp1(:) + cos(2.0*PI/3.0)*warp2(:) + cos(4.0*PI/3.0)*warp3(:)
  y(:) =  0.0*warp1(:) + sin(2.0*PI/3.0)*warp2(:) + sin(4.0*PI/3.0)*warp3(:)

  deallocate(xgll)
  deallocate(dl)
end subroutine shiftwarp

subroutine xyzToRst(x,y,z,r,s,t)
  ! INPUT: node coordinates in tetrahedron x,y,z
  ! DOES: transform from equilateral tetrahedron into reference tetrahedron
  ! OUTPUT: node coordinates in reference tetrahedron r,s,t
  implicit none
  real(kind=RKIND), dimension(:), intent(out) :: r,s,t
  real(kind=RKIND), dimension(:), intent(in) :: x,y,z
  real(kind=RKIND), dimension(3) :: v1,v2,v3,v4
  real(kind=RKIND), dimension(3,3) :: A
  real(kind=RKIND), dimension(3,Np) :: B,C

  A(:,:) = 0d0
  B(:,:) = 0d0
  C(:,:) = 0d0

  ! set vertices of the tetra
  v1 = (/-1d0, -1d0/dsqrt(3d0), -1d0/dsqrt(6d0) /)
  v2 = (/ 1d0, -1d0/dsqrt(3d0), -1d0/dsqrt(6d0) /)
  v3 = (/ 0d0,  2d0/dsqrt(3d0), -1d0/dsqrt(6d0) /)
  v4 = (/ 0d0,            0d0,   3d0/dsqrt(6d0) /)

  B(1,:) = x(:) - 0.5d0*(v2(1)+v3(1)+v4(1)-v1(1)) 
  B(2,:) = y(:) - 0.5d0*(v2(2)+v3(2)+v4(2)-v1(2)) 
  B(3,:) = z(:) - 0.5d0*(v2(3)+v3(3)+v4(3)-v1(3)) 

  A(:,1) = 0.5d0*(v2-v1)
  A(:,2) = 0.5d0*(v3-v1)
  A(:,3) = 0.5d0*(v4-v1)

  call invert(A)
  C = matmul(A,B)

  r= C(1,:)
  s= C(2,:)
  t= C(3,:)
end subroutine xyzToRst

subroutine rstToAbc(r,s,t,a,b,c)
  ! INPUT: node coordinates in reference tet r,s,t
  ! DOES: transform from reference tet into reference quad
  ! OUTPUT: node coordinates in reference quad a,b,c
  implicit none
  real(kind=RKIND), dimension(:), intent(in) :: r,s,t
  real(kind=RKIND), dimension(:), intent(out) :: a,b,c
  integer :: i

  do i=1,size(r)
    if (abs(s(i) + t(i)) > EPS) then
      a(i) = 2.0*(1.0+r(i))/(-s(i)-t(i))-1.0
    else
      a(i) = -1.0
    end if
    if (abs(t(i) - 1.) > EPS) then
      b(i) = 2.0*(1.0+s(i))/(1.0-t(i))-1.0
    else
      b(i) = -1.0
    end if
  end do
  c(:)=t(:)
end subroutine rstToAbc

subroutine rsToAb(r,s,a,b)
  ! transform from reference triange into reference quad
  implicit none
  real(kind=RKIND), dimension(:), intent(in) :: r,s
  real(kind=RKIND), dimension(:), intent(out) :: a,b
  integer :: i

  do i=1,size(r)
    if (abs(s(i) - 1.) > EPS) then
      a(i) = 2.0*(1.0+r(i))/(1.0-s(i))-1.0
    else
      a(i) = -1.0
    end if
  end do
  b(:)=s(:)
end subroutine rsToAb

subroutine xyToRs(x,y,r,s)
  ! transform from equilateral tri into reference tri
  implicit none
  real(kind=RKIND), dimension(:), intent(out) :: r,s
  real(kind=RKIND), dimension(:), intent(in) :: x,y
  real(kind=RKIND), dimension(Np) :: L1,L2,L3
  !
  L1(:) = (sqrt(3.0)*y(:)+1.0)/3.0
  L2(:) = (-3.0*x(:) - sqrt(3.0) * y(:) + 2.0)/6.0
  L3(:) = ( 3.0*x(:) - sqrt(3.0) * y(:) + 2.0)/6.0
  r(:) = -L2(:) + L3(:) - L1(:)
  s(:) = -L2(:) - L3(:) + L1(:)
end subroutine xyToRs

function getGll()
  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLL) :: x,GetGll
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL):: V
  real(kind=CUSTOM_REAL) :: xold
  integer :: i,k,p

  p = NGLL-1

  do i=1,NGLL
    x(i) = cos(pi*(i-1)/p)
  end do

! Newton-Raphson iteration
  do i=1,NGLL
    xold = 2.
    do while (abs(x(i)-xold) > eps)
      xold = x(i)
      V(i,1) = 1.
      V(i,2) = x(i)
      do k=2,p
        V(i,k+1) = ((2.*k-1.)*x(i)*V(i,k)-(k-1.)*V(i,k-1))/k
      end do
      x(i) = xold-(x(i)*V(i,NGLL)-V(i,p))/(NGLL*V(i,NGLL))
    end do
  end do
  GetGLL=x
end function

end module
