c      module simmetrix_snap
c      use iso_c_binding
c
c      interface
c        subroutine sim_closest_pos ( dx, dy, dz, id )
c     &    bind(C, NAME='sim_closest_pos')
c        use iso_c_binding
c          real(c_double), value :: dx, dy, dz
c          integer(c_int), value :: id
c        end subroutine
c
c        subroutine get_model_velocity (v)
c     &    bind(C, NAME='get_model_velocity')
c        use iso_c_binding
c          real(c_double) :: v
c        end subroutine
c      end interface
c
c      end module
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
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
c        use iso_c_binding
c        use simmetrix_snap
c
        include "common.h"
c
        real*8    x(numnp,nsd),  disp_snap(numnp,nsd)
        dimension iBC(nshg),     BC(nshg, 4)
        integer   surfID_snap(nshg)
        real*8    x_tmp(nsd),    x_crt(nsd)
        real*8    mag,           rad
        integer   i
c
        do i = 1, nshg
c... if surf ID is snapSurfID
          if (surfID_snap(i) .eq. snapSurfID) then
c... current x
            x_tmp(:) = x(i,:) + disp_snap(i,:)
c... correct x
c            call sim_closest_pos (disp_snap(i,1), disp_snap(i,2),
c     &                            disp_snap(i,3), i)
c... hardcoding; we should use simmetrix V_movedParamPoint to do this
            mag = sqrt(x_tmp(2)*x_tmp(2) + x_tmp(3)*x_tmp(3)) ! cylinder axis is x
            rad = 0.05
            x_crt(2) = x_tmp(2) / mag * rad
            x_crt(3) = x_tmp(3) / mag * rad
c... end hardcoding
c... update iBC and BC
            iBC(i) = ibset(iBC(i), 14)
            iBC(i) = ibset(iBC(i), 15)
            iBC(i) = ibset(iBC(i), 16)
            BC(i,1)= disp_snap(i,1)
            BC(i,2)= x_crt(2) - x(i,2)
            BC(i,3)= x_crt(3) - x(i,3)
          endif ! if equal to snapSurfID
        enddo
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine updateTDComp1BC ( x,       shpb,
     &                               iBC,     BC)
c
         use pointer_data
c
        include "common.h"
c
        real*8    x(numnp,nsd)
        dimension iBC(nshg),       BC(nshg,4),
     &            shpb(MAXTOP,maxsh,MAXQPT)
        integer   surfID_TDComp1(nshg)
        real*8    normal(nshg, nsd)
c
        if ( timeDepComp1Flag .eq. 1) then ! double check
c
        write(*,*) "Modify comp1_elas BC on surfID",timeDepComp1ID
c
        surfID_TDComp1 = zero
        normal = zero
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
          call calc_TDComp1(x,  shpb, mienb(iblk)%p,  miBCB(iblk)%p,
     &                      surfID_TDComp1, normal)
c
        enddo boundary_blocks ! end loop the boundary elements
c
        call resetTDComp1BC(x, surfID_TDComp1, normal, iBC, BC)
c
        endif ! end check TDComp1 flag
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine calc_TDComp1 (x,    shpb,     ienb,    iBCB,
     &                           surfID_TDComp1, normal)
c
        include "common.h"
c
        integer   surfID_TDComp1(nshg)
        dimension x(numnp,nsd),       shpb(nshl,ngaussb),
     &            iBCB(npro,ndiBCB),  ienb(npro,nshl)
        real*8    normal(nshg, nsd)
        integer   calc_factor(npro)
c
c... collect surf ID = timeDepComp1ID
c
        calc_factor(:) = 0
        do inode = 1, npro
          if (iBCB(inode,2) .eq. timeDepComp1ID) then
            calc_factor(inode) = 1
            do i = 1, nenbl ! only loop over vtx on the boundary
              surfID_TDComp1(ienb(inode,i)) = iBCB(inode,2)
            enddo
          endif
        enddo
c
c.... assemble the normal vector
c
        call calc_normal(x, shpb, calc_factor, ienb, iBCB, normal)
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine resetTDComp1BC (x, surfID_TDComp1, normal, iBC, BC)
c
        include "common.h"
c
        real*8    x(numnp,nsd)
        integer   surfID_TDComp1(nshg)
        real*8    normal(nshg, nsd)
        dimension iBC(nshg),     BC(nshg, 4)
        integer   maxDir(1)
c
        do i = 1, nshg
c... if surf ID is timeDepComp1ID
          if (surfID_TDComp1(i) .eq. timeDepComp1ID) then
c... if BC code is comp1 (iBC = 1 or 2 or 4)
            if ((ibits(iBC(i),14,3) .eq. 1) .or.
     &          (ibits(iBC(i),14,3) .eq. 2) .or.
     &          (ibits(iBC(i),14,3) .eq. 4)) then
c... calculate normal and normalize it
              maxDir = maxloc(abs(normal(i,:)))
              select case (maxDir(1))
c... if x of normal is the max
              case (1)
              iBC(i) = ibset(iBC(i), 14)
              iBC(i) = ibclr(iBC(i), 15)
              iBC(i) = ibclr(iBC(i), 16)
              BC(i,2)= normal(i,2) / normal(i,1)
              BC(i,3)= normal(i,3) / normal(i,1)
c... if y of normal is the max
              case (2)
              iBC(i) = ibclr(iBC(i), 14)
              iBC(i) = ibset(iBC(i), 15)
              iBC(i) = ibclr(iBC(i), 16)
              BC(i,2)= normal(i,1) / normal(i,2)
              BC(i,3)= normal(i,3) / normal(i,2)
c... if z of normal is the max
              case (3)
              iBC(i) = ibclr(iBC(i), 14)
              iBC(i) = ibclr(iBC(i), 15)
              iBC(i) = ibset(iBC(i), 16)
              BC(i,2)= normal(i,1) / normal(i,3)
              BC(i,3)= normal(i,2) / normal(i,3)
              end select
            endif ! if BC code is comp1
          endif ! if equal to timeDepComp1ID
        enddo
c
        return
        end
c
