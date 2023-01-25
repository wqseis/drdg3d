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
!along with NEXD 2D. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------

module mod_exchange
  use mpi
  use mod_para,  only : RKIND, CUSTOM_REAL, &
                        Nfp, dimens
  use mod_mpi,   only : isendV_real,        &
                        irecV_real,         &
                        wait_req
  use mod_types, only : meshvar, buffvar

  implicit none

  contains

  subroutine init_buff(mesh,buff)
    implicit none

    type(meshvar) :: mesh
    type(buffvar) :: buff

    if(mesh%nproc < 2) return

    allocate(buff%q_send(Nfp*mesh%mpi_ne*dimens,mesh%mpi_nn))
    allocate(buff%q_rec(Nfp*mesh%mpi_ne*dimens,mesh%mpi_nn))
    allocate(buff%qi(Nfp,dimens,mesh%mpi_ne,mesh%mpi_nn))

  end subroutine

  subroutine exchange_data(mesh,u,buff)
    implicit none
    type(meshvar) :: mesh
    type(buffvar) :: buff
    real(kind=RKIND),allocatable,dimension(:,:) :: u

    integer :: i,j,k,ie,ierr,tag,dest,req,req_r,c,je,face
    !integer :: myrank,nproc,sizes
    integer :: sizes

    if(mesh%nproc < 2) return

    ! mpi comunication
    ! build send buffer
    !qm = q
    tag = 100
    buff%q_send = 0.
    do i=1,mesh%mpi_nn
      do ie=1,mesh%mpi_ne ! loop over interface elements
        je   = mesh%mpi_connection(i,ie,1)
        face = mesh%mpi_connection(i,ie,2)
        if (je > 0) then
          do k=1,dimens
            !tmpu = reshape(u(:,je,k),(/NGLL,NGLL/))
            !call get_face1(tmpu,face,tmpuface)
            !tmpuface = u(mesh%vmapM(:,face),je,k)
            do j = 1,Nfp
              buff%q_send((ie-1)*dimens*Nfp + (k-1)*Nfp + j,i) = &
                u(mesh%vmapM(j,face,je),k)
            end do
          end do
        end if
      end do
    end do ! all interfaces

    !if (.false.) then
    !send and rec
    do i=1,mesh%mpi_nn
      dest=mesh%mpi_neighbor(i)-1
      sizes = (mesh%mpi_ne*dimens*Nfp)
      call isendV_real(buff%q_send(:,i),sizes,dest,tag,req,CUSTOM_REAL)
      call irecV_real(buff%q_rec(:,i),sizes,dest,tag,req_r,CUSTOM_REAL)
      call wait_req(req)
      call wait_req(req_r)
    end do

    !endif

    !call MPI_Barrier(MPI_COMM_WORLD,ierr)

    !unpack mpi buffers
    do i=1,mesh%mpi_nn
      c = 1
      do ie=1,mesh%mpi_ne
        do k=1,dimens
          do j=1,Nfp
            if ( mesh%mpi_connection(i,ie,1) > 0) then
              buff%qi(j,k,ie,i) = buff%q_rec(c,i)
            end if
            c = c+1
          end do
        end do
      end do
    end do

    !call sync_mpi()
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine

end module
