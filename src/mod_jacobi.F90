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

module mod_jacobi

  use mod_para,   only : RKIND, NGLL, Np
  use mod_matrix, only : invert
  use mod_nodes,  only : rstoab, rsttoabc

  implicit none

contains

subroutine jacobiP(P,x,alpha,beta,N)
  ! INPUT: x,,alpha,beta,N: Parameters of Jacobi polynomials
  ! DOES: Calculates the Jacobi polynomials
  ! OUTPUT: P: Array with Jacobi polynomial
  implicit none
  real(kind=RKIND), dimension(:), intent(in) :: x
  integer, intent(in) :: N
  real(kind=RKIND), dimension(N+1,size(x)) :: PL
  real(kind=RKIND), dimension(:), intent(out) :: P
  real(kind=RKIND) :: gamma0, gamma1, aold, anew, h1, bnew
  real(kind=RKIND), intent(in) :: alpha, beta
  integer :: dim, i

  dim=size(x)
  gamma0 = 2**(alpha+beta+1.0)/(alpha+beta+1.0)* &
  !lacz_gamma(dble(alpha+1.0))*lacz_gamma(dble(beta+1.0))/lacz_gamma(dble(alpha+beta+1.0))
  !PL(1,:) = 1d0/dsqrt(gamma0)
  lacz_gamma((alpha+1.0))*lacz_gamma((beta+1.0))/lacz_gamma((alpha+beta+1.0))
  PL(1,:) = 1d0/sqrt(gamma0)

  if (N == 0) then
    P=PL(1,:)
    return
  end if

  gamma1 = (alpha+1.0)*(beta+1.0)/(alpha+beta+3.0)*gamma0
  !PL(2,:) = ((alpha+beta+2.0)*x(:)/2d0 + (alpha-beta)/2d0) / dsqrt(gamma1)
  PL(2,:) = ((alpha+beta+2.0)*x(:)/2d0 + (alpha-beta)/2d0) / sqrt(gamma1)

  if (N == 1) then
    P=PL(N+1,:)
    return
  end if

  ! repeat value in recurrence
  aold = 2d0/(2d0+alpha+beta)*dsqrt((alpha+1d0)*(beta+1d0)/(alpha+beta+3d0))

  ! forward recurrence using the symetry of the recurrence
  do i=1,N-1
    h1 = 2d0*i+alpha+beta
    anew = 2d0/(h1+2d0)*dsqrt( (i+1d0) * (i+1d0+alpha+beta)*(i+1d0+alpha)*(i+1d0+beta)/(h1+1d0)/(h1+3d0) )
    bnew = - (alpha**2-beta**2)/h1/(h1+2d0)
    PL(i+2,:) = 1d0/anew*( -aold*PL(i,:) + (x(:)-bnew)*PL(i+1,:))
    aold=anew
  end do

  P=PL(N+1,:)
end subroutine jacobiP

subroutine gradJacobiP(dP,x,alpha,beta,N)
  ! INPUT: x,,alpha,beta,N: Parameters of Jacobi gradient
  ! DOES: compute grad of jacobi polynomials
  ! OUTPUT: P: Array with Jacobi gradient
  implicit none
  real(kind=RKIND), dimension(:), intent(out) :: dP
  real(kind=RKIND), dimension(:), intent(in) :: x
  real(kind=RKIND), intent(in) :: alpha, beta
  real(kind=RKIND), dimension(size(dp)) ::h1
  integer, intent(in) :: N
  if (N==0) then
    dP(:)=0
  else
    call jacobiP(h1,x,alpha+1,beta+1,N-1)
    dp(:) = dsqrt(N*(N+alpha+beta+1d0))*h1(:)
  end if
end subroutine gradJacobiP

recursive function lacz_gamma(a) result(g)
  implicit none
  real(kind=RKIND), intent(in) :: a
  real(kind=RKIND) :: g

  real(kind=RKIND), parameter :: pi = 3.141592653589793238463d0
  integer, parameter :: cg = 7

  ! these precomputed values are taken by the sample code in Wikipedia,
  ! and the sample itself takes them from the GNU Scientific Library
  real(kind=RKIND), dimension(0:8), parameter :: p = &
       (/ 0.99999999999980993d0, 676.5203681218851d0, -1259.1392167224028d0, &
       771.32342877765313d0, -176.61502916214059d0, 12.507343278686905d0, &
       -0.13857109526572012d0, 9.9843695780195716d-6, 1.5056327351493116d-7 /)

  real(kind=RKIND) :: t, w, x
  integer :: i

  x = a

  if ( x < 0.5d0 ) then
    !g = pi / ( sin(pi*x) * lacz_gamma(1d0-x) )
    g = pi / ( sin(pi*x) * lacz_gamma(1.0-x) )
  else
    x = x - 1d0
    t = p(0)
    do i=1, cg+1
      t = t + p(i)/(x+real(i))
    end do
    w = x + dble(cg) + 0.5d0
    !g = dsqrt(2d0*pi) * w**(x+0.5d0) * dexp(-w) * t
    g = dsqrt(2d0*pi) * w**(x+0.5d0) * exp(-w) * t
  end if
end function lacz_gamma

subroutine Vdm2D(v,r,s)
  implicit none
  ! initialize the 2D vandermonde matrix
  real(kind=RKIND), dimension(:), intent(in) :: r,s
  real(kind=RKIND), dimension(size(r),size(r)), intent(out) :: v
  real(kind=RKIND), dimension(size(r)) :: a,b
  integer :: i,j,sk
  !
  call rstoab(r,s,a,b)
  sk=1
  do i=0,NGLL-1
    do j=0,NGLL-1-i
      call simplex2DP(v(:,sk),a,b,i,j)
      sk=sk+1
    end do
  end do
end subroutine Vdm2D

subroutine invVdm2D(v,w,trans,n)
  ! calulates v^-1' if trans==0 -> inverse if trans == 1 => inverse transponierte
  integer, intent(in) ::n
  real(kind=RKIND), dimension(:,:), intent(in) :: v
  real(kind=RKIND), dimension(size(v(:,1)),n) :: w
  integer, intent(in) :: trans

  w=v
  call invert(w)

  if (trans==1) then
    w=transpose(w)
  end if

end subroutine invVdm2D

subroutine Vdm3D(v,r,s,t)
  ! INPUT: node coordinates in reference tetrahedron r,s,t
  ! DOES: initialize the 3D vandermonde matrix
  ! OUTPUT: van-der-monde matrix v
  implicit none
  real(kind=RKIND), dimension(:), intent(in) :: r,s,t
  real(kind=RKIND), dimension(size(r),Np), intent(out) :: v
  real(kind=RKIND), dimension(size(r)) :: a,b,c
  integer :: i,j,k,sk
  !
  call rstToAbc(r,s,t,a,b,c)  ! transform into reference tet

  sk=1
  do i=0,NGLL-1
    do j=0,NGLL-1-i
      do k=0,NGLL-1-i-j
        call simplex3DP(v(:,sk),a,b,c,i,j,k)
        sk=sk+1
      end do
    end do
  end do
end subroutine Vdm3D

subroutine invVdm3D(v,w,trans)
  ! INPUT: van-der-monde matrix v, flag trans
  ! DOES: calulates v^-1', if trans==0 -> inverse; if trans == 1 => inverse transpose
  ! OUTPUT: inverse (transpose) van-der-monde matrix w
  real(kind=RKIND), dimension(:,:), intent(in) :: v
  real(kind=RKIND), dimension(size(v(:,1)),Np) :: w
  integer, intent(in) :: trans

  w=v
  call invert(w)
  if (trans==1) then
    w=transpose(w)
  end if
end subroutine invVdm3D

subroutine gradVdm3D(v3dr,v3ds,v3dt,r,s,t)
  ! INPUT: node coordinates reference tet r,s,t
  ! DOES: build the gradient of the modal basis (i,j,k) at (r,s,t) in referenze tet
  ! OUTPUT: gradients to basis r,s,t (v3dr, v3ds, v3dt)
  implicit none
  real(kind=RKIND), dimension(:), intent(in) :: r,s,t
  real(kind=RKIND), dimension(size(r),Np), intent(out) :: v3dr,v3ds,v3dt
  real(kind=RKIND), dimension(size(r)) :: a,b,c
  integer :: i,j,k,sk
  !
  call rstToabc(r,s,t,a,b,c)

  sk=1
  do i=0,NGLL-1
    do j=0,NGLL-1-i
      do k=0,NGLL-1-i-j
        call gradSimplex3DP(v3dr(:,sk),v3ds(:,sk),v3dt(:,sk),a,b,c,i,j,k)
        sk=sk+1
      end do
    end do
  end do
end subroutine gradVdm3D

subroutine simplex2DP(Pr,a,b,i,j)
  implicit none
  real(kind=RKIND), dimension(:), intent(out) :: Pr
  real(kind=RKIND), dimension(:), intent(in):: a,b
  real(kind=RKIND), dimension(size(a)) :: h1,h2
  integer, intent(in) :: i,j
  real(kind=RKIND) :: zero,idvar1

  zero=0d0
  idvar1=2d0*i+1d0

  call jacobiP(h1,a,zero,zero,i)
  call jacobiP(h2,b,idvar1,zero,j)

  Pr(:) = dsqrt(2d0)*h1*h2*(1d0-b)**i
end subroutine simplex2DP

subroutine simplex3DP(Pr,a,b,c,i,j,k)
  ! INPUT: node coordinates in reference quad a,b,c, orders i,j,k
  ! DOES: Calculate orthonormal polynomial on simplex at a,b,c
  ! RETURNS: orthonormal polynomial Pr
  real(kind=RKIND), dimension(:), intent(out) :: Pr
  real(kind=RKIND), dimension(:), intent(in):: a,b,c
  real(kind=RKIND), dimension(size(a)) :: h1,h2,h3
  integer, intent(in) :: i,j,k
  real(kind=RKIND) :: zero,idvar1,idvar2

  zero=0d0
  idvar1=2d0*i+1d0
  idvar2=2d0*(i+j)+2d0

  call jacobiP(h1,a,zero,zero,i)
  call jacobiP(h2,b,idvar1,zero,j)
  call jacobiP(h3,c,idvar2,zero,k)

  Pr(:) = 2d0*dsqrt(2d0)*h1*h2*((1d0-b)**i)*h3*((1d0-c)**(i+j))
end subroutine simplex3DP

subroutine gradSimplex2DP(dPdr,dPds,a,b,id,jd)
  ! computes the differentiation matrices Dr and Ds
  implicit none
  real(kind=RKIND), dimension(:), intent(in) ::a,b
  integer, intent(in) :: id,jd
  real(kind=RKIND), dimension(size(a)), intent(out) :: dPdr,dPds
  real(kind=RKIND), dimension(size(a)) :: fa,dfa
  real(kind=RKIND), dimension(size(b)) :: gb,dgb
  real(kind=RKIND), dimension(size(a)) :: temp
  real(kind=RKIND) :: zero,idvar1

  zero=0d0
  idvar1=2d0*id+1d0

  call jacobiP(fa,a,zero,zero,id)
  call gradJacobiP(dfa,a,zero,zero,id)
  call jacobiP(gb,b,idvar1,zero,jd)
  call gradJacobiP(dgb,b,idvar1,zero,jd)

  dPdr(:)=dfa(:)*gb(:)
  if (id>0) then
    dPdr(:)=dPdr(:)*((0.5d0*(1d0-b(:)))**(id-1d0))
  end if

  dPds(:) = dfa(:)*(gb(:)*(0.5d0*(1d0+a(:))))
  if (id>0) then
    dPds(:)=dPds(:)*((0.5d0*(1d0-b(:)))**(id-1d0))
  end if

  temp(:) = dgb(:)*((0.5d0*(1d0-b(:)))**id)
  if (id>0) then
    temp(:)=temp(:)-0.5d0*id*gb(:)*((0.5d0*(1d0-b(:)))**(id-1d0))
  end if
  dPds(:) = dPds(:)+fa(:)*temp(:)

  !normalize
  dPdr(:) = 2d0**(id+0.5d0)*dPdr(:)
  dPds(:) = 2d0**(id+0.5d0)*dPds(:)
end subroutine gradSimplex2DP

subroutine gradSimplex3DP(dPdr,dPds,dPdt,a,b,c,id,jd,kd)
  ! computes the differentiation matrices Dr and Ds
  implicit none
  real(kind=RKIND), dimension(:), intent(in) ::a,b,c
  integer, intent(in) :: id,jd,kd
  real(kind=RKIND), dimension(size(a)), intent(out) :: dPdr,dPds,dPdt
  real(kind=RKIND), dimension(size(a)) :: fa,dfa
  real(kind=RKIND), dimension(size(b)) :: gb,dgb
  real(kind=RKIND), dimension(size(c)) :: hc,dhc
  real(kind=RKIND), dimension(size(a)) :: temp
  real(kind=RKIND) :: zero,idvar1,idvar2

  zero=0d0
  idvar1=2d0*id+1d0
  idvar2=2d0*(id+jd)+2d0

  call jacobiP(fa,a,zero,zero,id)
  call gradJacobiP(dfa,a,zero,zero,id)
  call jacobiP(gb,b,idvar1,zero,jd)
  call gradJacobiP(dgb,b,idvar1,zero,jd)
  call jacobiP(hc,c,idvar2,zero,kd)
  call gradJacobiP(dhc,c,idvar2,zero,kd)

  ! r-deriverate

  dPdr(:) = dfa(:)*(gb(:)*hc(:))
  if (id>0) then
    dPdr(:) = dPdr(:)*((0.5d0*(1d0-b(:)))**(id-1d0))
  end if
  if ((id+jd)>0) then
    dPdr(:) = dPdr(:)*((0.5d0*(1d0-c(:)))**(id+jd-1d0))
  end if

  ! s-deriverate

  dPds(:) = 0.5d0*(1d0+a(:))*dPdr(:)
  temp(:) = dgb(:)*((0.5d0*(1d0-b(:)))**id)
  if (id>0) then
    temp(:) = temp(:)+(-0.5d0*id)*(gb(:)*(0.5d0*(1d0-b(:)))**(id-1d0))
  end if
  if ((id+jd)>0) then
    temp(:) = temp(:)*((0.5d0*(1d0-c(:)))**(id+jd-1d0))
  end if
  temp(:) = fa(:)*temp(:)*hc(:)
  dPds(:) = dpds(:)+temp(:)

  ! t-deriverate

  dPdt(:) = 0.5d0 * (1d0+a(:)) * dPdr(:) + 0.5d0*(1d0+b(:))*temp(:)
  temp(:) = dhc(:) *((0.5d0*(1d0-c(:)))**(id+jd))
  if ((id+jd)>0) then
    temp(:) = temp(:) -0.5d0*(id+jd)*(hc*((0.5d0*(1d0-c(:)))**(id+jd-1d0)))
  end if
  temp(:) = fa(:) * (gb(:)*temp(:))
  temp(:) = temp(:) *((0.5d0*(1d0-b(:)))**id)

  dPdt(:) = dPdt(:) + temp(:)

  !normalize
  dPdr(:) = (2d0**(2d0*id+jd+1.5d0))*dPdr(:)
  dPds(:) = (2d0**(2d0*id+jd+1.5d0))*dPds(:)
  dPdt(:) = (2d0**(2d0*id+jd+1.5d0))*dPdt(:)
end subroutine gradSimplex3DP

subroutine Dmatrices3D(dr,ds,dt,r,s,t,v)
  ! INPUT: node coordinates reference tet r,s,t, van-der-monde-matrix v
  ! DOES: calculates the differentiation matrices on the simplex at (r,s,t)
  ! OUTPUT: differentiation matrices dr,ds,dt with respect to r,s,t
  implicit none
  real(kind=RKIND), dimension(:), intent(in) :: r,s,t
  real(kind=RKIND), dimension(:,:), intent(in) :: v
  real(kind=RKIND), dimension(size(v(:,1)),Np) :: dr,ds,dt,w,vr,vs,vt

  call gradVdm3D(vr,vs,vt,r,s,t)
  call invVdm3D(v,w,0)

  dr=matmul(vr,w)
  ds=matmul(vs,w)
  dt=matmul(vt,w)

end subroutine Dmatrices3D

end module
