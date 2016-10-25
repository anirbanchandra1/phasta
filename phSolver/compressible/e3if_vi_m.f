      module e3if_vi_m
c
        use e3if_inp_m
        use e3if_param_m
c
        implicit none
c
      contains
c
      subroutine calc_vi(p,n,u)
c
        real*8, pointer, intent(in) :: p(:), n(:,:), u(:,:)
c
        integer :: isd
        real*8 :: vimag,v1,t1,c1
c
        select case (vi_model)
        case (const_vi)
          vimag = vi_mag
        case default
          call error ('ERROR in e3if_vi:',' vi_model is not supported.',vi_model)
        end select
c
        select case (vi_ramping)
        case (no_ramp)
          c1 = vimag
        case (linear_ramp)
          t1 = ramp_time
          v1 = vimag
c
          if (time <= t1) then
            c1 = v1*time/t1
          else
            c1 = v1
          endif
        case default
          call error ('ERROR in e3if_vi:',' vi_ramping is not supported,',vi_ramping)
        end select
c
        do isd = 1,nsd
          vi(:,isd) = c1*n(:,isd) + u(:,isd)
        enddo
c
      end subroutine calc_vi
c
      subroutine calc_vi_area_node(sum_vi_area_l,shp,WdetJif,nshl)
c
        real*8, dimension(:,:,:), intent(inout) :: sum_vi_area_l
        real*8, dimension(:,:),   intent(in)    :: shp
        real*8, dimension(:),     intent(in)    :: WdetJif
        integer, intent(in) :: nshl
c
        integer :: n,isd
c
        do n = 1,nshl
          do isd = 1,nsd
            sum_vi_area_l(:,n,isd) = sum_vi_area_l(:,n,isd) 
     &                             + shp(:,n)*vi(:,isd)*area(:)
c     &                             + shp(:,n)*vi(:,isd)*WdetJif(:)
          enddo
          sum_vi_area_l(:,n,nsd+1) = sum_vi_area_l(:,n,nsd+1) + shp(:,n)*area(:)
c          sum_vi_area_l(:,n,nsd+1) = sum_vi_area_l(:,n,nsd+1) + shp(:,n)*WdetJif(:)
        enddo
c
      end subroutine calc_vi_area_node
c
      end module e3if_vi_m
