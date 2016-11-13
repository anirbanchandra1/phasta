      subroutine calc_kappa_error(x,iel0,nelblkif,nsd,nshg)
c
#define case 1 
c
        use pointer_data
        use workfc_m
        use if_global_m
c
        implicit none
c
        include "mpif.h"
        real*8, intent(in) :: x(nshg,nsd)
        integer, intent(in) :: iel0(nelblkif)
        integer, intent(in) :: nsd,nshg,nelblkif
c
        integer :: inode(6), n, iel, iblk, npro, nsum, my_nsum, ierr
        integer, dimension(:,:), pointer :: ienif0, ienif1
        real*8 :: this_kappa, my_err, err, this_err, kappa_exact, max_err, my_max_err
        logical, dimension(nshg) :: check
c
        my_err = 0.d0
        my_max_err = 0.d0
        my_nsum = 0
        check = .false.
c
        do iblk = 1,nelblkif
c
          npro = iel0(iblk+1) - iel0(iblk)
          ienif0 => mienif0(iblk)%p
          ienif1 => mienif1(iblk)%p
c
          do iel = 1,npro
            inode(1:3) = ienif0(iel,1:3)
            inode(4:6) = ienif1(iel,1:3)
            do n = 1,6
              if (.not. check(inode(n))) then
                check(inode(n)) = .true.
              else
                cycle
              endif
#if case == 1
              kappa_exact = 10.d0
#endif
              my_nsum = my_nsum + 1
              this_err = (norm2(if_kappa(inode(n),1:nsd)) - kappa_exact)**2
              my_err = my_err + this_err
              my_max_err = max(my_max_err,this_err)
c
            enddo
          enddo
        enddo
c
        call mpi_allreduce(my_err,err,1,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
        call mpi_allreduce(my_max_err,max_err,1,MPI_DOUBLE_PRECISION,MPI_MAX, MPI_COMM_WORLD,ierr)
        call mpi_allreduce(my_nsum,nsum,1,MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD,ierr)
c
        err = sqrt(err) / real(nsum,8)
c
        if (myrank == 0) then
          write(*,100) nsum,err,max_err
        endif
c
100   format(' Total interface nodes, L2, L_inf norms: ',i6,x,2(2x,e24.16))
c
      end subroutine calc_kappa_error