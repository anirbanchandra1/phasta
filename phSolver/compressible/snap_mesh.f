        subroutine updateSnapSurfBC ( x,       disp_snap,
     &                                iBC,     BC)
c
         use pointer_data
c
        include "common.h"
c
        real*8    x(numnp,nsd),    disp_snap(numnp,nsd)
        dimension iBC(nshg),       BC(nshg,4)
        integer   surfID_snap(nshg)
c
        if ( snapSurfFlag .eq. 1) then ! double check
c
        surfID_snap = zero
c
c.... loop over the boundary elements
c
        boundary_blocks: do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel
c
          call checkSnapSurfID(mienb(iblk)%p, miBCB(iblk)%p, surfID_snap)
c
        enddo boundary_blocks ! end loop the boundary elements
c
        call resetSnapBC(x, disp_snap, iBC, BC, surfID_snap)
c
        endif ! end check snap flag
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine checkSnapSurfID (ienb, iBCB, surfID_snap)
c
        include "common.h"
c
        integer   surfID_snap(nshg)
        dimension iBCB(npro,ndiBCB),  ienb(npro,nshl)
c
        do inode = 1, npro
          if (iBCB(inode,2) .eq. snapSurfID) then
            do i = 1, nenbl ! only loop over vtx on the boundary
              surfID_snap(ienb(inode,i)) = iBCB(inode,2)
            enddo
          endif
        enddo
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine resetSnapBC ( x,       disp_snap,
     &                           iBC,     BC,     surfID_snap)
c
        include "common.h"
c
        real*8    x(numnp,nsd),  disp_snap(numnp,nsd)
        dimension iBC(nshg),     BC(nshg, 4)
        integer   surfID_snap(nshg)
        real*8    x_tmp(nsd),    x_crt(nsd)
        real*8    mag,           rad
c
        do i = 1, nshg
c... if surf ID is snapSurfID
          if (surfID_snap(i) .eq. snapSurfID) then
c... current x
            x_tmp(:) = x(i,:) + disp_snap(i,:)
c... correct x
c... hardcoding; we should use simmetrix V_movedParamPoint to do this
            mag = sqrt(x_tmp(2)*x_tmp(2) + x_tmp(3)*x_tmp(3)) ! cylinder axis is x
            rad = 0.05
            x_crt(1) = x_tmp(1)
            x_crt(2) = x_tmp(2) / mag * rad
            x_crt(3) = x_tmp(3) / mag * rad
c... end hardcoding
c... update iBC and BC
            iBC(i) = ibset(iBC(i), 14)
            iBC(i) = ibset(iBC(i), 15)
            iBC(i) = ibset(iBC(i), 16)
            BC(i,1)= x_crt(1) - x(i,1)
            BC(i,2)= x_crt(2) - x(i,2)
            BC(i,3)= x_crt(3) - x(i,3)
          endif ! if equal to snapSurfID
        enddo
c
        return
        end


