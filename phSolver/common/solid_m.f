      module solid_m
c
        use pointer_data
c
        implicit none
c
        integer :: iblk_solid
c
        type (r3d), dimension(MAXBLK2) ::  b !for solid,left Cauchy_green tensor,added
        type (r3d), dimension(MAXBLK2) ::  b_dot!time derivative of b,added
        type (r3d), dimension(MAXBLK2) ::  b_af! b at time step n+af,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b !left Cauchy_green tensor on the boundary,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b_dot!time derivative of b on the boudary,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b_af! b at time step n+af on the boudary,added
c
        integer ::iblk_now
c
      contains
c
        subroutine alloc_solid
c
        end subroutine alloc_solid
c
        subroutine free_solid
c
        end subroutine free_solid
c
      end module solid_m
