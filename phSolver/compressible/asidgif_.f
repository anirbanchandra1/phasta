      subroutine asidgif_
     & (
     &   nshg_, nshl0_, nshl1_, nenl0_, nenl1_, lcsyst0_, lcsyst1_,
     &   nflow_, npro_, ndof_, nsd_, ipord_, numnp_, nqpt_,
     &   gbytes_, sbytes_, flops_,
     &   res,
     &   egmassif00, egmassif01, egmassif10, egmassif11,
     &   y,        x,       umesh,
     &   shpif0,   shpif1,  shgif0,  shgif1,
     &   qwtif0,   qwtif1,
     &   ienif0,   ienif1,
     &   materif0, materif1,
     &   time_,
     &   sum_vi_area
     & )  
         use hierarchic_m
         use local_m
         use e3if_m
c
            implicit none
c
            integer, intent(in) :: nshg_, nshl0_, nshl1_, nenl0_, nenl1_,lcsyst0_,lcsyst1_
            integer, intent(in) :: nflow_, npro_, ndof_, nsd_, ipord_, numnp_, nqpt_
            integer, intent(inout) :: gbytes_, sbytes_, flops_
            real*8, dimension(nshg_,nflow_), intent(inout) :: res
            real*8, dimension(:,:,:), pointer, intent(inout) :: egmassif00,egmassif01,egmassif10,egmassif11
            real*8, dimension(nshg_,ndof_),  intent(in)    :: y
            real*8, dimension(nshg_,nsd_),   intent(in)    :: x
            real*8, dimension(numnp_, nsd_), intent(inout) :: umesh
            real*8, dimension(nshl0_,nqpt_),intent(in)   :: shpif0
            real*8, dimension(nshl1_,nqpt_),intent(in)   :: shpif1
            real*8, dimension(nsd_,nshl0_,nqpt_), intent(in)  :: shgif0
            real*8, dimension(nsd_,nshl1_,nqpt_), intent(in)  :: shgif1
            real*8, dimension(nqpt_), intent(in) :: qwtif0, qwtif1
            integer, dimension(:,:), pointer, intent(in)   :: ienif0, ienif1
            integer, intent(in)   :: materif0, materif1
            real*8, intent(in) :: time_
            real*8, intent(inout) :: sum_vi_area(:,:)
c
      integer :: i, iel, i0,i1,n
      real*8 :: r0,r1,n0(3),n1(3),t_ramp,c2
      real*8, parameter :: rmin = 0.01d0, rmax = 0.9d0,
     &                     vint = -1.0d0
c
            call setparam_e3if
     &      (
     &        nshg_,nshl0_,nshl1_,nenl0_,nenl1_,lcsyst0_,lcsyst1_,
     &        npro_,ndof_,nsd_,nflow_,ipord_,nqpt_,
     &        egmassif00,egmassif01,egmassif10,egmassif11,
     &        materif0, materif1,
     &        time_
     &      )
c
          call malloc_e3if
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
          if (ipord_ .gt. 1) then
           call getsgn(ienif0,sgn0,nshl0,nenl0)
           call getsgn(ienif1,sgn1,nshl1,nenl1)
        endif
c
c... localize
c
        call localy(y, ycl0, ienif0, ndof, 'gather  ', nshg, nshl0, npro, ipord, gbytes_)
        call localy(y, ycl1, ienif1, ndof, 'gather  ', nshg, nshl1, npro, ipord, gbytes_)
c
        call localx(x, xl0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro, gbytes_)
        call localx(x, xl1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro, gbytes_)
c
        call localx(umesh, umeshl0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro, gbytes_)
        call localx(umesh, umeshl1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro, gbytes_)
c
       call e3if(shpif0,shpif1,shgif0,shgif1,qwtif0,qwtif1)
c
c.... assemble the local residual arrays
c
c      do iel = 1,npro
c        do n = 1,nshl0
c          i0 = ienif0(iel,n)
c          write(*,200) myrank,0,iel,n,i0,x(i0,:)
c        enddo
c        do n = 1,nshl1
c          i1 = ienif1(iel,n)
c          write(*,200) myrank,1,iel,n,i1,x(i1,:)
c        enddo
c      enddo
c
c      write(*,*) 'CHECK IN ASIDGIF BEFORE LOCAL:'
c      do iel = 1,npro
c        do n = 1,nshl0
c          i0 = ienif0(iel,n)
c          write(*,100) myrank,0,iel,n,i0,rl0(iel,n,:)
c          write(*,100) myrank,0,iel,n,i0,res(i0,:)
c        enddo
c        do n = 1,nshl1
c          i1 = ienif1(iel,n)
c          write(*,100) myrank,1,iel,n,i1,rl1(iel,n,:)
c          write(*,100) myrank,1,iel,n,i1,res(i1,:)
c        enddo
c      enddo
        call local (res, rl0, ienif0, nflow, 'scatter ', nshg,nshl0,npro,ipord,sbytes_,flops_)
        call local (res, rl1, ienif1, nflow, 'scatter ', nshg,nshl1,npro,ipord,sbytes_,flops_)
c      write(*,*) 'CHECK IN ASIDGIF AFTER LOCAL:'
      do iel = 1,npro
        do n = 1,nshl0
          i0 = ienif0(iel,n)
c          write(*,100) myrank,0,iel,n,i0,res(i0,:)
        enddo
        do n = 1,nshl1
          i1 = ienif1(iel,n)
c          write(*,100) myrank,1,iel,n,i1,res(i1,:)
        enddo
      enddo
100   format('[',i2,'] iel,n,ienif',i1,': ',i4,i2,i6,5e24.16)
200   format('[',i2,'] iel,n,ienif',i1,': ',i4,i2,i6,3e24.16)
c
        call local (sum_vi_area, sum_vi_area_l0, ienif0, nsd+1, 'scatter ', nshg, nshl0,npro,ipord,sbytes_,flops_)
        call local (sum_vi_area, sum_vi_area_l1, ienif1, nsd+1, 'scatter ', nshg, nshl1,npro,ipord,sbytes_,flops_)
c
        call free_e3if
c
      end subroutine asidgif_
