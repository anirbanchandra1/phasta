c-----------------------------------------------------------------------
c
c    Predict solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrPredictElas (disp,umesh,dt)
      
      include "common.h"
      
      real*8   disp(nshg,nelas),umesh(numnp,nsd),dt
c
c      disp = zero
      disp = umesh*dt
c
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrectElas (xold, x, disp, Dy)
      
      include "common.h"
      
      real*8   xold(numnp,nsd), x(numnp,nsd)
      real*8   disp(nshg,nelas), Dy(nshg,nelas)
c      
      disp = disp + Dy
      x    = xold + disp
c
      return
      end

      subroutine itrUpdateElas (xold, x)
      
      include "common.h"
      
      real*8   xold(numnp,nsd), x(numnp,nsd)
c      
      xold    = x
c
      return
      end

