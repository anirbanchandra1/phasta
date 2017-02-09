      module genshp_m
c
        use global_const_m
        use number_def_m
        use elmpar_m
        use shpdat_m
        use genpar_m
        use blkdat_m
        use intpt_m
c
        implicit none
c
        contains

        subroutine genshpif
     &(            shpif, shpif0, shpif1
     &,            shgif, shgif0, shgif1
     &)
c
          real*8, dimension(MAXTOP,      maxsh, MAXQPT), intent(out) :: shpif   ! shape function
          real*8, dimension(MAXTOP,      maxsh, MAXQPT), intent(out) :: shpif0  ! shape function          for interface element0
          real*8, dimension(MAXTOP,      maxsh, MAXQPT), intent(out) :: shpif1  ! shape function          for interface element1
          real*8, dimension(MAXTOP, nsd, maxsh, MAXQPT), intent(out) :: shgif   ! shape function gradient 
          real*8, dimension(MAXTOP, nsd, maxsh, MAXQPT), intent(out) :: shgif0  ! shape function gradient 
          real*8, dimension(MAXTOP, nsd, maxsh, MAXQPT), intent(out) :: shgif1  ! shape function gradient 
c
          integer :: i, iblk, id, itpif
c
c... loop over the interface blocks and 
c    generate shape functions and gradients
c
          do iblk = 1, nelblif
c
            nshl0   = lcblkif(iblkif_nshl0,iblk)
            nshl1   = lcblkif(iblkif_nshl1,iblk)
            id      = lcblkif(iblkif_topology,iblk)
c
            do i = 1, nintif(id)
              call shpTet(ipord, Qptif(itp_tet,1:3,i), shpif(itp_tet,:,i), shgif(itp_tet,:,:,i))
              call shp6w (ipord, Qptif(itp_wedge,1:3,i), shpif(itp_wedge,:,i),shgif(itp_wedge,:,:,i))
            enddo
c
            shgif(itp_tet,:,:,1:nintif(id)) = 
     &        shgif(itp_tet,:,:,1:nintif(id)) / two 
c
            cycle
c NOTE: need to remove below!!!
c....
c
            select case (id)
            case (itpif_tet_tet)
c
              do i = 1, nintif(id)
                call shpTet(ipord, Qptif(id,1:3,i), shpif0(itp_tet,:,i), shgif0(itp_tet,:,:,i))
                call shpTet(ipord, Qptif(id,1:3,i), shpif1(itp_tet,:,i), shgif1(itp_tet,:,:,i))
              enddo
c
              shgif0(itp_tet,:,1:nshl0,1:nintif(id)) = 
     &          shgif0(itp_tet,:,1:nshl0,1:nintif(id)) / two 
              shgif1(itp_tet,:,1:nshl1,1:nintif(id)) = 
     &          shgif1(itp_tet,:,1:nshl1,1:nintif(id)) / two 
c
            case (itpif_tet_wedge)
c
              do i = 1, nintif(id)
                call shpTet(ipord,Qptif(id,1:3,i),shpif0(itp_tet,:,i),shgif0(itp_tet,:,:,i))
                call shp6w (ipord,Qptif(id,1:3,i),shpif1(itp_wedge,:,i),shgif1(itp_wedge,:,:,i))
              enddo
c
              shgif0(itp_tet,:,1:nshl0,1:nintif(id)) = 
     &          shgif0(itp_tet,:,1:nshl0,1:nintif(id)) / two 
c
            case (itpif_wedge_tet)
c
              do i = 1, nintif(id)
                call shp6w (ipord,Qptif(id,1:3,i),shpif0(itp_wedge,:,i),shgif0(itp_wedge,:,:,i))
                call shpTet(ipord,Qptif(id,1:3,i),shpif1(itp_tet,:,i),shgif1(itp_tet,:,:,i))
              enddo
c
              shgif1(itp_tet,:,1:nshl1,1:nintif(id)) = 
     &          shgif1(itp_tet,:,1:nshl1,1:nintif(id)) / two 
c
            case (itpif_wedge_wedge)
c
              do i = 1, nintif(id)
                call shp6w (ipord,Qptif(id,1:3,i),shpif0(itp_wedge,:,i),shgif0(itp_wedge,:,:,i))
                call shp6w (ipord,Qptif(id,1:3,i),shpif1(itp_wedge,:,i),shgif1(itp_wedge,:,:,i))
              enddo
c
            case default
c
c.... nonexistent element
c
              call error ('genshpif  ', 'elem Cat', id)
c
            end select
c
          enddo
c
        end subroutine genshpif
c
      end module genshp_m
