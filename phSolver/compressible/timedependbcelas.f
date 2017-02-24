       subroutine timeDependBCElas(x, iBC, BC, BC_flow, umeshold)
c
c-----------------------------------------------------------------
c
c-----------------------------------------------------------------
c
      include "common.h"
      include "mpif.h"
c
c
      real*8    x(numnp,nsd)
      real*8    disp(numnp,nsd)
      real*8    umeshold(numnp,nsd)
      dimension iBC(nshg),        BC(nshg,3), BC_flow(nshg,3)
      integer   casenumber
      real*8    acc
      real*8    totalForce(3),    objMass
      dimension Forin(3),         Forout(3)
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
c.... bcast force to all processors
        if (numpe > 1) then
          do j = 1, nsrfCM
            Forin  = (/ Force(1,j), Force(2,j), Force(3,j) /)
            call MPI_BCAST (Forin(1), 3, MPI_DOUBLE_PRECISION,
     &                      master,      MPI_COMM_WORLD,ierr)
            Force(1:3,j) = Forin(1:3)
          enddo
        endif
c.... collect total force
        do j = 1,nsrfCM
          totalForce(:) = totalForce(:) + Force(:,j)
        enddo ! end collect total force

        acc = totalForce(1) / objMass

        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (abs(x(i,1)-0.8) .lt. 0.79) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 0.004225)) then
c.... 0.004225 = 0.065^2
            BC(i,1)   = umeshold(i,1) + acc * Delt(1)
            BC(i,2)   = zero
            BC(i,3)   = zero

c.... >>> hard-coding update flow velocity boundry condition on mortar surfaces
            BC_flow(i,1:3) = BC(i,1:3)
c.... <<< hard-coding

          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 2
c
c.... Update BC value based on geom and iBC
c.... rotate + trans + shrink
c
      if ( casenumber .eq. 3 ) then
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (abs(x(i,1)) .lt. 0.2) .and.
     &         (abs(x(i,2)) .lt. 0.2) .and.
     &         (abs(x(i,3)) .lt. 0.2) ) then
c
            disp(i,1) = -0.005 * (x(i,1)) !+ 0.02
     &                + (x(i,1)) * (cos(pi/600) - 1.0)
     &                - (x(i,2)) *  sin(pi/600)
            disp(i,2) = !-0.005 * (x(i,2)) !+ 0.01
     &                + (x(i,2)) * (cos(pi/600) - 1.0)
     &                + (x(i,1)) *  sin(pi/600)
            disp(i,3) = -0.005 * x(i,3)
c
            BC(i,:)   = disp(i,:) / Delt(1)
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 3
c
      return
      end
c

