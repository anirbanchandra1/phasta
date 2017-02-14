      subroutine shpMeasure (x, ien, meshq, tetVolm, errorcount)
c---------------------------------------------------------
c This subroutine calculate the shape quality 
c based on the mean ratio method. 
c
c---------------------------------------------------------
c
        include "common.h"
c
        real*8    meshq(npro) 
        dimension x(numnp, nsd),        ien(npro, nshl),
     &            xl(npro, nenl, nsd)
c
        dimension tetEdge(npro, 6, nsd), 
     &            tetVolm(npro), 
     &            tetSumEdge(npro)
c
        integer errorcount(2)
c
c---------------------------------------------------------
c          
       call localx(x,      xl,     ien,    nsd,    'gather  ')
c
       tetEdge(:,1,:) = xl(:,1,:) - xl(:,4,:)
       tetEdge(:,2,:) = xl(:,2,:) - xl(:,4,:)
       tetEdge(:,3,:) = xl(:,3,:) - xl(:,4,:)
       tetEdge(:,4,:) = xl(:,1,:) - xl(:,3,:)
       tetEdge(:,5,:) = xl(:,2,:) - xl(:,3,:)
       tetEdge(:,6,:) = xl(:,1,:) - xl(:,2,:)
c
       tetVolm(:) = 1.0/6.0 * ( 
     &       tetEdge(:,1,1) * ( tetEdge(:,2,2) * tetEdge(:,3,3)
     &                        - tetEdge(:,2,3) * tetEdge(:,3,2) )
     &     - tetEdge(:,1,2) * ( tetEdge(:,2,1) * tetEdge(:,3,3)
     &                        - tetEdge(:,2,3) * tetEdge(:,3,1) )
     &     + tetEdge(:,1,3) * ( tetEdge(:,2,1) * tetEdge(:,3,2)
     &                        - tetEdge(:,2,2) * tetEdge(:,3,1) ) )
c
       tetSumEdge(:) = tetEdge(:,1,1) * tetEdge(:,1,1)
     &               + tetEdge(:,1,2) * tetEdge(:,1,2)
     &               + tetEdge(:,1,3) * tetEdge(:,1,3)
     &               + tetEdge(:,2,1) * tetEdge(:,2,1)
     &               + tetEdge(:,2,2) * tetEdge(:,2,2)
     &               + tetEdge(:,2,3) * tetEdge(:,2,3)
     &               + tetEdge(:,3,1) * tetEdge(:,3,1)
     &               + tetEdge(:,3,2) * tetEdge(:,3,2)
     &               + tetEdge(:,3,3) * tetEdge(:,3,3)
     &               + tetEdge(:,4,1) * tetEdge(:,4,1)
     &               + tetEdge(:,4,2) * tetEdge(:,4,2)
     &               + tetEdge(:,4,3) * tetEdge(:,4,3)
     &               + tetEdge(:,5,1) * tetEdge(:,5,1)
     &               + tetEdge(:,5,2) * tetEdge(:,5,2)
     &               + tetEdge(:,5,3) * tetEdge(:,5,3)
     &               + tetEdge(:,6,1) * tetEdge(:,6,1)
     &               + tetEdge(:,6,2) * tetEdge(:,6,2)
     &               + tetEdge(:,6,3) * tetEdge(:,6,3)
c
       meshq(:) = 15552.0 * tetVolm(:) * abs(tetVolm(:))
     &          / ( tetSumEdge(:) * tetSumEdge(:) * tetSumEdge(:) ) 
c
       do i = 1, npro
          if (meshq(i) .gt. 1.0000001 .and. printflag1 .eq. 0) then
             errorcount(1) = errorcount(1) + 1
          endif
          if (meshq(i) .lt. 0.0 .and. printflag2 .eq. 0) then 
             tetVolm(i) = abs (tetVolm(i))
             errorcount(2) = errorcount(2) + 1
          endif
       enddo
c
c.... return
c
       return
       end
