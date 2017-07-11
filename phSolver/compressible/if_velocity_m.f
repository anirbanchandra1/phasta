      module if_velocity_m
c
c----------------------------------------
c    aims to calculate the interface
c    velocity at the global level
c----------------------------------------
        use workfc_m
        use number_def_m
        use pointer_data
        use blkdat_m
c
        implicit none
c
        real*8, pointer :: sum_vi_area(:,:)    ! interface velocity weighted by interface area
c
      contains
c
      subroutine init_sum_vi_area(nshg,nsd)
        integer, intent(in) :: nshg,nsd
        if (associated(sum_vi_area)) 
     &    deallocate (sum_vi_area)
        allocate (sum_vi_area(nshg,nsd+1))
        sum_vi_area(:,1:nsd) = zero
        sum_vi_area(:,1+nsd) = one
      end subroutine init_sum_vi_area
c
      subroutine destruct_sum_vi_area
        if (associated(sum_vi_area))
     &    deallocate(sum_vi_area)
      end subroutine destruct_sum_vi_area
c
      subroutine set_if_velocity 
     & (
     &  BC, iBC, umesh, disp, x, dt, ilwork,
     &  nshg, ndofBC, nsd, nelblif, nlwork
     & )
c
        include "mpif.h"
c
        real*8,  intent(inout) ::  BC(nshg,3)
        integer, intent(inout) :: iBC(nshg)
        real*8,  dimension(nshg,nsd), intent(inout)    :: umesh, disp
        real*8,  dimension(nshg,nsd), intent(in)    :: x
        real*8, intent(in) :: dt
        integer, intent(in)    :: ilwork(nlwork)
        integer, intent(in) :: nshg, ndofBC, nsd, nelblif, nlwork
c
        integer :: iblk, iel, npro,inode, i0, i1, n, ierr
        integer, pointer :: ienif0(:,:), ienif1(:,:)

        if (numpe > 1) then
          call commu (sum_vi_area(:,1:3), ilwork, nsd, 'in ')
          call commu (sum_vi_area(:,4), ilwork, 1, 'in ')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        endif
c 
        if (numpe > 1) then
          call commu (sum_vi_area(:,1:3), ilwork, nsd, 'out')
          call commu (sum_vi_area(:,4), ilwork, 1, 'out')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        endif
c
c        do inode = 1,nshg
c
c ... NOT SURE IF THIS IS THE BEST IF :
c
c          if (abs(sum_vi_area(inode,nsd+1)) > zero) then
c            umesh(inode,:) = sum_vi_area(inode,:) / sum_vi_area(inode,nsd+1)
c            BC(inode,:) = umesh(inode,:)
c      write(*,100) 'AFTER: ', myrank,inode, x(inode,:), sum_vi_area(inode,:),umesh(inode,:)
c      write(*,200) 'AFTER: ', myrank,inode, umesh(inode,:)
c          endif
c        enddo
c
	!call itrBCvi(sum_vi_area(:,1:3),ac, iBC, BC, iper, ilwork) 
	open (unit = 2000, file = "debug.dat")
	!HACK-AC!
	sum_vi_area(106,2)=0
	sum_vi_area(106,3)=0
	sum_vi_area(107,2)=0
	sum_vi_area(107,3)=0
	sum_vi_area(108,2)=0
        sum_vi_area(108,3)=0
        sum_vi_area(109,2)=0
        sum_vi_area(109,3)=0
	sum_vi_area(5,2)=0
        sum_vi_area(5,3)=0
        sum_vi_area(6,2)=0
        sum_vi_area(6,3)=0
	sum_vi_area(7,2)=0
        sum_vi_area(7,3)=0
        sum_vi_area(8,2)=0
        sum_vi_area(8,3)=0

        do iblk = 1,nelblif
          npro  = lcblkif(1,iblk+1) - lcblkif(1,iblk)
          ienif0 => mienif0(iblk)%p
          ienif1 => mienif1(iblk)%p
          do iel = 1,npro
            do n = 1,3  ! only triangles on the surface
              i0 = ienif0(iel,n)
              i1 = ienif1(iel,n)
              umesh(i0,:) = sum_vi_area(i0,:) / sum_vi_area(i0,nsd+1)
              umesh(i1,:) = sum_vi_area(i1,:) / sum_vi_area(i1,nsd+1)
c              BC(i0,:) = sum_vi_area(i0,:) / sum_vi_area(i0,nsd+1)
c              BC(i1,:) = sum_vi_area(i1,:) / sum_vi_area(i1,nsd+1)
              BC(i0,:) = umesh(i0,:) * dt
              BC(i1,:) = umesh(i1,:) * dt
      write(2000,200) 'AFTER: ', myrank,i0, umesh(i0,:)
      write(2000,200) 'AFTER: ', myrank,i1, umesh(i1,:)
            enddo
          enddo
        enddo
c
100   format(a,'[',i2,'] ',i6,3f7.3,x,7e14.6)
200   format(a,'[',i2,'] ',i6,3e14.6)
      end subroutine set_if_velocity
c
      end module if_velocity_m
