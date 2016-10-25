      module solid_data_m
c
        use pointer_data
        implicit none
c
        type (r3d), dimension(MAXBLK2) ::  b !for solid,left Cauchy_green tensor,added
        type (r3d), dimension(MAXBLK2) ::  b_dot!time derivative of b,added
        type (r3d), dimension(MAXBLK2) ::  b_af! b at time step n+af,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b !left Cauchy_green tensor on the boundary,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b_dot!time derivative of b on the boudary,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b_af! b at time step n+af on the boudary,added
c
        integer, pointer :: is_solid(:)
        integer, parameter :: b_size = 6
c
        type solid_t
          logical :: is_active
          integer :: nel
          real*8, pointer :: b(:)
        end type solid_t
c
        type(solid_t) :: solid_p
c
      end module solid_data_m
c
      module solid_m
c
        use solid_data_m
        implicit none
c
      contains
c
        subroutine alloc_solid
          use matdat_def_m
          use number_def_m
          use elmpar_m
          use blkdat_m
          use intpt_m
          implicit none
c
          integer :: mattype, iblk, npro
c
c.... allocate space for solid arrays 
 
          blocks_loop: do iblk = 1, nelblk
c
            mattype = lcblk(i_mattype,iblk)
            lcsyst = lcblk(3,iblk)
            ngauss = nint(lcsyst)
            npro = lcblk(1,iblk+1) - lcblk(1,iblk)
c
            if (mat_eos(mattype,1).eq.ieos_solid_1)then
c
              allocate (b(iblk)%p(npro,ngauss,b_size))
              allocate (b_dot(iblk)%p(npro,ngauss,b_size))
              allocate (b_af(iblk)%p(npro,ngauss,b_size))
c
              if (associated(solid_p%b)) then
c
              else
c
                b(iblk)%p(:,:,:)= one
                b_af(iblk)%p(:,:,:) = one
c
              endif
            else
              cycle   
            endif
c
          enddo blocks_loop
c
          boundary_blocks_loop: do iblk = 1, nelblb
c
            mattype = lcblkb(i_mattype,iblk)
            lcsyst = lcblkb(3,iblk)
            ngauss = nintb(lcsyst)
            npro = lcblkb(1,iblk+1) - lcblkb(1,iblk)
c
            if (mat_eos(mattype,1).eq.ieos_solid_1)then
c
              allocate (bdy_b(iblk)%p(npro,ngaussb,b_size))
              allocate (bdy_b_dot(iblk)%p(npro,ngaussb,b_size))
              allocate (bdy_b_af(iblk)%p(npro,ngaussb,b_size))
c
              if (associated(solid_p%b)) then
c
              else
c
                bdy_b(iblk)%p(:,:,:)= one
                bdy_b_af(iblk)%p(:,:,:) = one
c
              endif
            else
              cycle
            endif
c..
          enddo boundary_blocks_loop
c
        end subroutine alloc_solid
c
        subroutine free_solid
c
        end subroutine free_solid
c
      end module solid_m
