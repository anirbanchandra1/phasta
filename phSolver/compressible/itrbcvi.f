        subroutine itrBCvi (actual_vi, iBC, BC)
c--apply flow BCs on sum_vi_area
c
c----------------------------------------------------------------------
c
c This program satisfies the flow BCs on the actual_vi variables
c
c input:
c  y      (nshg,nflow)   : y variables 
c  iBC    (nshg)        : Boundary Condition Code
c  BC     (nshg,ndofBC) : boundary condition constraint parameters
c  ylimit (3,nflow)     : (1,:) limiting flag
c                         (2,:) lower bound
c                         (3,:) upper bound
c output:
c  y      (nshg,nflow)   : Adjusted V value(s) corresponding to a 
c                           constraint d.o.f.
c  umesh  (numnp,nsd)    : mesh velocity. FOR ALE 
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension actual_vi(nshg,3),             iBC(nshg),
     &            BC(nshg,ndofBC)
          
c.... ------------------------->  Velocity  <--------------------------
c.... 3D
c
c.... x1-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 1)
            actual_vi(:,1) =  BC(:,3)  - BC(:,4) * actual_vi(:,2)
     &                         - BC(:,5) * actual_vi(:,3)
          endwhere
c
c.... x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 2)
            actual_vi(:,2) = BC(:,3)  - BC(:,4) * actual_vi(:,1)
     &                        - BC(:,5) * actual_vi(:,3)
          endwhere
c
c.... x1-velocity and x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 3)
            actual_vi(:,1) =  BC(:,3)  - BC(:,4) * actual_vi(:,3)
            actual_vi(:,2) =  BC(:,5)  - BC(:,6) * actual_vi(:,3)
          endwhere
c
c.... x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 4)
            actual_vi(:,3) = BC(:,3) - BC(:,4) * actual_vi(:,1)
     &                       - BC(:,5) * actual_vi(:,2)
          endwhere
c
c.... x1-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 5)
            actual_vi(:,1) = BC(:,3) - BC(:,4) * actual_vi(:,2)
            actual_vi(:,3) = BC(:,5) - BC(:,6) * actual_vi(:,2)
          endwhere
c
c.... x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 6)
            actual_vi(:,2) = BC(:,3)  - BC(:,4) * actual_vi(:,1)
            actual_vi(:,3) = BC(:,5)  - BC(:,6) * actual_vi(:,1)
          endwhere
c
c.... x1-velocity, x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 7)
            actual_vi(:,1) =  BC(:,3)
            actual_vi(:,2) =  BC(:,4)
            actual_vi(:,3) =  BC(:,5) 
          endwhere
c
c       endif
c
c.... end of velocity
c
        return
        end subroutine itrBCvi

