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

module mod_source

  implicit none

contains

function svf_ricker(t,fc,tdelay)

  implicit none
  real*8 :: svf_ricker
  real*8 :: f0,r,rr,t,fc,tdelay
  real*8,parameter :: pi = 3.141593

  f0 = sqrt(pi)/2.
  r = pi * fc * (t-tdelay)
  rr = r**2

  svf_ricker = r*(3.-2.*rr)*exp(-rr)*f0*pi*fc

end function


end module
