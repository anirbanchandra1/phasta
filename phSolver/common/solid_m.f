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
      contains
c
        subroutine alloc_solid(lcblk, lcblkb, nelblk, nelblb,
     &                         nint,  nintb,  MAXTOP, i_iniSolid) !check
        use matdat_def_m
        use number_def_m
        implicit none
c
        integer, dimension (10,nelblk+1) :: lcblk
        integer, dimension (10,nelblb+1) :: lcblkb
        integer, dimension (MAXTOP) :: nint, nintb
        integer :: nelblk, nelblb, MAXTOP, i_iniSolid
c
        integer :: mattyp_s, lcsyst_s, npro_s, ngauss_s
        integer :: mattyp_sb, lcsyst_sb, npro_sb, ngauss_sb
        integer :: iblk
c
c.... allocate space for solid arrays 
 
       blocks_loop: do iblk = 1, nelblk
         mattyp_s = lcblk(7,iblk)
         lcsyst_s = lcblk(3,iblk)
         ngauss_s = nint(lcsyst_s)
         npro_s = lcblk(1,iblk+1) - lcblk(1,iblk)
!for solid block only
         if (mat_eos(mattyp_s,1).eq.ieos_solid_1)then
           allocate (b(iblk)%p(npro_s,ngauss_s,6))!for solid
           allocate (b_dot(iblk)%p(npro_s,ngauss_s,6)) !for solid,added
           allocate (b_af(iblk)%p(npro_s,ngauss_s,6)) !for solid,added
c......... these arrays need initialization
           if (i_iniSolid .eq. 1)then
              b(iblk)%p    = zero 
              b_dot(iblk)%p = zero
              b(iblk)%p(:,:,1)= one
              b(iblk)%p(:,:,2)= one
              b(iblk)%p(:,:,3)= one
              b_af(iblk)%p = b(iblk)%p
           endif
c
          endif
c..
      enddo blocks_loop
c
c
c
c 
      boundary_ blocks_loop: do iblk = 1, nelblb
         mattyp_sb = lcblkb(7,iblk)
         lcsyst_sb = lcblkb(3,iblk)
         ngauss_sb = nintb(lcsyst_sb)
         npro_sb = lcblkb(1,iblk+1) - lcblkb(1,iblk)
!for solid block only
         if (mat_eos(mattyp_sb,1).eq.ieos_solid_1)then
           allocate (bdy_b(iblk)%p(npro_sb,ngauss_sb,6))!for solid
           allocate (bdy_b_dot(iblk)%p(npro_sb,ngauss_sb,6)) !for solid,added
           allocate (bdy_b_af(iblk)%p(npro_sb,ngauss_sb,6)) !for solid,added
c......... these arrays need initialization
           if (i_iniSolid .eq. 1)then
             bdy_b(iblk)%p    = zero
             bdy_b_dot(iblk)%p = zero
             bdy_b(iblk)%p(:,:,1)= one
             bdy_b(iblk)%p(:,:,2)= one
             bdy_b(iblk)%p(:,:,3)= one
             bdy_b_af(iblk)%p = bdy_b(iblk)%p
           endif
c
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
