      subroutine asidgif
     & (
     &   res,
     &   y,        x,       umesh,
     &   shpif0,   shpif1,  shgif0,  shgif1,
     &   qwtif0,   qwtif1,
     &   ienif0,   ienif1,
     &   sum_vi_area, if_normals
     & )  
          use hierarchic_m
          use local_m
          use e3if_m
          use e3if_geom_m
          use if_global_m
          use conpar_m
c
          implicit none
c
          real*8, dimension(nshg,nflow), intent(inout) :: res
          real*8, dimension(nshg,ndof),  intent(in)    :: y
          real*8, dimension(nshg,nsd),   intent(in)    :: x
          real*8, dimension(nshg, nsd), intent(inout) :: umesh
          real*8, dimension(nshl0,nqpt),intent(in)   :: shpif0
          real*8, dimension(nshl1,nqpt),intent(in)   :: shpif1
          real*8, dimension(nsd,nshl0,nqpt), intent(in)  :: shgif0
          real*8, dimension(nsd,nshl1,nqpt), intent(in)  :: shgif1
          real*8, dimension(nqpt), intent(in) :: qwtif0, qwtif1
          integer, dimension(:,:), pointer, intent(in)   :: ienif0, ienif1
          real*8, pointer, intent(inout) :: sum_vi_area(:,:)
          real*8, pointer, intent(in) :: if_normals(:,:)
      integer :: i0,i1,iel,n,i
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
          if (ipord .gt. 1) then
           call getsgn(ienif0,sgn0,nshl0,nenl0)
           call getsgn(ienif1,sgn1,nshl1,nenl1)
        endif
c
c... localize
c
        call localy(y, ycl0, ienif0, ndof, 'gather  ', nshg, nshl0, npro, ipord)
        call localy(y, ycl1, ienif1, ndof, 'gather  ', nshg, nshl1, npro, ipord)
c
        call localx(x, xl0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
        call localx(x, xl1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
c
        call localx(umesh, umeshl0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
        call localx(umesh, umeshl1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
c
        call localx(if_normals, if_normal_l0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
        call localx(if_normals, if_normal_l1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
c
        if (associated(if_kappa)) then
          call localx(if_kappa, if_kappa_l0,  ienif0, nsd, 'gather  ', nshg, nenl0, npro)
          call localx(if_kappa, if_kappa_l1,  ienif1, nsd, 'gather  ', nshg, nenl1, npro)
        endif
c
        call e3if(shpif0,shpif1,shgif0,shgif1,qwtif0,qwtif1)
c
c.... assemble the local residual arrays
c
      i0 = ienif0(1,1)
      i1 = ienif1(2,1)
      !write(*,*) 'i0 before: ',i0,res(i0,:)
      !write(*,10) 'i1 before: ',i1,res(i1,:)
        call local (res, rl0, ienif0, nflow, 'scatter ', nshg,nshl0,npro,ipord)
        call local (res, rl1, ienif1, nflow, 'scatter ', nshg,nshl1,npro,ipord)
      !write(*,*) 'i0 after: ',i0,res(i0,:)
      do iel = 1,npro
        do n = 1,nshl1
          i0 = ienif0(iel,n)
          i1 = ienif1(iel,n)
           !write(*,10) 'rl0: ',n,ienif0(iel,n),rl0(iel,n,:)
           !write(*,10) 'rl1: ',n,ienif1(iel,n),rl1(iel,n,:)
           !write(*,10) 'i0 after: ',n,i0,res(i0,:)
           !write(*,10) 'i1 after: ',n,i1,res(i1,:)
        enddo
      enddo
      !write(*,10) 'i1 after: ',i1,res(i1,:)
10    format(a8,2i6,5e24.16)
c
        call local (sum_vi_area, sum_vi_area_l0, ienif0, nsd+1, 'scatter ', nshg, nshl0,npro,ipord)
        call local (sum_vi_area, sum_vi_area_l1, ienif1, nsd+1, 'scatter ', nshg, nshl1,npro,ipord)
c
      end subroutine asidgif
