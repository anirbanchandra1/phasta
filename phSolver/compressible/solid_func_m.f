      module solid_func_m
c
        use e3_def_m
        use solid_m
        use eqn_state_m
c
        implicit none
c
        real*8, dimension(:,:), pointer :: dudx, dudy, dudz
        real*8, dimension(:,:,:), pointer :: AS
c
      contains
c
      subroutine calc_solid
c
        allocate(AS(npro,6,6))
c
        call calc_as_matrix
c        call setB_af
c        call get_det
c        call get_det
c
        deallocate(AS)
c
      end subroutine calc_solid
c
      subroutine calc_as_matrix
c
      end subroutine calc_as_matrix
c
      end module solid_func_m
