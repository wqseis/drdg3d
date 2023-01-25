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

module mod_solve

  use mod_para, only : RKIND

  implicit none

  integer, parameter :: sp = selected_real_kind(9,49) !< Short precision
  integer, parameter :: lp = selected_real_kind(15,99) !< Long precision
  !integer, parameter :: wp = lp !< Working precision
  integer, parameter :: wp = RKIND !< Working precision

contains

subroutine Regula_Falsi(V,Phi,eta,sigma_n,psi,V0,a)
  ! a nonlinear solver for slip-rate using Regula-Falsi
  ! solve: V + sigma_n/eta*f(V,psi) - Phi/eta = 0

  implicit none
  !character(256),intent(in)::problem
  real(kind = wp), intent(inout) :: V          ! magnitude of slip-velocity (slip-rate)
  real(kind = wp), intent(in) :: Phi           ! stress transfer functions
  real(kind = wp), intent(in) :: eta           ! half of harmonic average of shear impedance
  real(kind = wp), intent(in) :: sigma_n       ! compressive normal stress
  real(kind = wp), intent(in) :: psi           ! state variable
  real(kind = wp), intent(in) :: a             ! direct effect parameter
  real(kind = wp), intent(in) :: V0            ! reference velocity
  real(kind = wp) :: fv, fl, fr, Vl, Vr, V1    ! function values
  real(kind = wp) :: err, tol
  integer ::  k, maxit

  err = 1.0_wp
  tol = 1.0e-12_wp

  Vl = 0.0_wp             ! lower bound
  Vr = Phi/eta            ! upper bound

  k = 1
  maxit = 5000

  call f(fv,V,Phi,eta,sigma_n,psi,V0,a)
  call f(fl,Vl,Phi,eta,sigma_n,psi,V0,a)
  call f(fr,Vr,Phi,eta,sigma_n,psi,V0,a)

  if (abs(fv) .le. tol) then ! check if the guess is a solution
    V = V
    return
  end if

  if (abs(fl) .le. tol) then ! check if the lower bound is a solution
    V = Vl
    return
  end if

  if (abs(fr) .le. tol) then  ! check if the upper bound is a solution
    V = Vr
    return
  end if

  ! else solve for V iteratively using regula-falsi
  do while((k .le. maxit) .and. (err .ge. tol))

    V1 = V

    if (fv*fl > tol) then
      Vl = V
      V = Vr - (Vr-Vl)*fr/(fr-fl)
      fl = fv
    elseif (fv*fr > tol) then
      Vr = V
      V = Vr -(Vr-Vl)*fr/(fr-fl)
      fr = fv
    end if

    call f(fv,V,Phi,eta,sigma_n,psi,V0,a)

    err = abs(V-V1)
    k = k + 1

  end do

  if (k .ge. maxit) print*, err, k, 'maximum iteration reached. solution didnt convergence!'

end subroutine Regula_Falsi

subroutine f(fv,V,Phi,eta,sigma_n,psi,V0,a)
  !                                                                   
  ! This subroutine evaluates the function:
  ! f:= V + sigma_n/eta*a*f(V,psi) - Phi/eta = 0
  ! where V = sqrt(v1**2 + v2**2)

  implicit none
  !character(256),intent(in)::problem
  real(kind = wp), intent(in) :: V             ! magnitude of slip-velocity (slip-rate)
  real(kind = wp), intent(in) :: Phi           ! stress transfer functions
  real(kind = wp), intent(in) :: eta           ! half of harmonic average of shear impedance
  real(kind = wp), intent(in) :: sigma_n       ! compressive normal stress
  real(kind = wp), intent(in) :: psi           ! state variable
  real(kind = wp), intent(in) :: a             ! direct effect parameter
  real(kind = wp), intent(in) :: V0            ! reference velocity
  real(kind = wp), intent(inout) ::fv          ! function values

  !select case(problem)

  !case('rate-weakening')

  fv = V + a*sigma_n/eta*asinh(0.5_wp*V/V0*exp(psi/a)) - Phi/eta

  !end select
end subroutine f


end module
