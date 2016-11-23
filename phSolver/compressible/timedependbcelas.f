       subroutine timeDependBCElas(x, iBC, BC, umeshold)
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
      include "common.h"
c
c
      real*8    x(numnp,nsd)
      real*8    disp(numnp,nsd)
      real*8    umeshold(numnp,nsd)
      dimension iBC(nshg),        BC(nshg,3)
      integer   casenumber
      real*8    acc
      real*8    totalForce(3),    objMass
c
      if (elasFDC .gt. 0) then
        casenumber = elasFDC
      else
        write(*,*) "Please use Force-driven as Mesh Elas Model"
      endif
c
c.... Update BC value based on geom and iBC
c
      if ( casenumber .eq. 1 ) then
        acc = 60000.0
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (abs(x(i,1)-0.8) .lt. 0.79) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 0.004225)) then
c.... 0.004225 = 0.065^2
            BC(i,1)   = acc * lstep * Delt(1)
            BC(i,2)   = zero
            BC(i,3)   = zero
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 1
c
c.... Update BC value based on total force on the object
c
      if ( casenumber .eq. 2 ) then
        totalForce(:) = zero
        objMass = 15.0 ! kg
        do j = 1,nsrfCM
          totalForce(:) = totalForce(:) + Force(:,j)
        enddo ! end collect total force

        acc = totalForce(1) / objMass
c.... debugging
        write(*,*) "projectile acc: ", acc
c.... end debugging

        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (abs(x(i,1)-0.8) .lt. 0.79) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 0.004225)) then
c.... 0.004225 = 0.065^2
            BC(i,1)   = umeshold(i,1) + acc * Delt(1)
            BC(i,2)   = zero
            BC(i,3)   = zero
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 2
c
      return
      end
c

