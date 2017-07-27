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
c.... dynamic origin, x translation, rotation frequence
      real*8    dyn_org,   xtsl,  rotf
c
      if (elasFDC .gt. 0) then
        write(*,*) "use Force-driven case:", elasFDC
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
            BC(i,1)   = acc * lstep * Delt(1) * Delt(1)
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
            BC(i,1)   = umeshold(i,1) + acc * Delt(1) * Delt(1)
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
c          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
          if ( (ibits(iBC(i),3,3) .eq. 7) .and.
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
            BC(i,:)   = disp(i,:)
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 3
c
c.... Update BC value based on geom and iBC
c.... twist one end of a cylinder
c
      if ( casenumber .eq. 4 ) then
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 0.501) .and.
     &         (x(i,1) .gt. 0.499) ) then
c
            disp(i,1) = 0.0
            disp(i,2) = x(i,2) * (cos(pi/100) - 1.0)
     &                - x(i,3) *  sin(pi/100)
            disp(i,3) = x(i,3) * (cos(pi/100) - 1.0)
     &                + x(i,2) *  sin(pi/100)
c
            BC(i,:)   = disp(i,:)
          endif ! end if inside box
        enddo ! end loop numnp
      endif  ! end case 4
c
c.... test case 5
c.... original mov+rot+shrk
c
      if ( casenumber .eq. 5 ) then
        if (lstep .eq. 0) then
        dyn_org    = 2.625e-3
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 6.65e-3) .and. (x(i,1) .gt. -1.4e-3) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 2.89e-6) ) then
            if( x(i,2) .ge. 0.6e-3 ) then ! top
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! top tail
                disp(i,1) = -0.05 * (x(i,1)-0.5e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.05 * x(i,3)
              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! top head
                disp(i,1) = -0.05 * (x(i,1)-4.75e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.05 * x(i,3)
              else ! top middle
                disp(i,1) = -0.05 * (x(i,1)-2.625e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
                disp(i,3) = -0.05 * x(i,3)
              endif ! end top--switch head middle tail
            else if( x(i,2) .le. -0.6e-3 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! bottom tail
                disp(i,1) = -0.05 * (x(i,1)-0.5e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.05 * x(i,3)
              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! bottom head
                disp(i,1) = -0.05 * (x(i,1)-4.75e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.05 * x(i,3)
              else ! bottom middle
                disp(i,1) = -0.05 * (x(i,1)-2.625e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
                disp(i,3) = -0.05 * x(i,3)
              endif ! end bottom--switch head middle tail
            else ! middle
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! middle tail
                disp(i,1) = -0.05 * (x(i,1)-0.5e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
                disp(i,3) = -0.05 * x(i,3)
              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! middle head
                disp(i,1) = -0.05 * (x(i,1)-4.75e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
                disp(i,3) = -0.05 * x(i,3)
              else ! middle middle
                disp(i,1) = -0.05 * (x(i,1)-2.625e-3) !+ 0.02
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
                disp(i,3) = -0.05 * x(i,3)
              endif ! end middle--switch head middle tail
            endif ! end if switch top middle bottom
            BC(i,1:3)   = disp(i,1:3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
        endif ! first time step
      endif ! end if case 5
c
c.... end test case 5
c
c
c.... test case 6
c.... original mov+rot+shrk
c
      if ( casenumber .eq. 6 ) then
        xtsl     = 1.0000000000000e-4
        dyn_org  = 2.6250000000000e-3 + lstep * xtsl
        write(*,*) "current lstep:", lstep, "dyn_org:", dyn_org
        rotf     = 180.0
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 6.65e-3) .and. (x(i,1) .gt. -1.4e-3) .and.
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 2.89e-6) ) then
            if( x(i,2) .ge. 0.6e-3 ) then ! top
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! top tail
                disp(i,1) = -0.05 * (x(i,1)-dyn_org+2.125e-3) + xtsl
     &                    + (x(i,1)-dyn_org+2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org+2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! top head
                disp(i,1) = -0.05 * (x(i,1)-dyn_org-2.125e-3) + xtsl
     &                    + (x(i,1)-dyn_org-2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org-2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              else ! top middle
                disp(i,1) = -0.05 * (x(i,1)-dyn_org) + xtsl
     &                    + (x(i,1)-dyn_org) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,2)-1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)-1.125e-3) !+ 0.01
     &                    + (x(i,2)-1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,1)-dyn_org) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -0.6e-3 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! bottom tail
                disp(i,1) = -0.05 * (x(i,1)-dyn_org+2.125e-3) + xtsl
     &                    - (x(i,1)-dyn_org+2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org+2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! bottom head
                disp(i,1) = -0.05 * (x(i,1)-dyn_org-2.125e-3) + xtsl
     &                    - (x(i,1)-dyn_org-2.125e-3) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org-2.125e-3) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              else ! bottom middle
                disp(i,1) = -0.05 * (x(i,1)-dyn_org) + xtsl
     &                    - (x(i,1)-dyn_org) * (cos(pi/rotf) - 1.0)
     &                    + (x(i,2)+1.125e-3) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)+1.125e-3) !- 0.01
     &                    - (x(i,2)+1.125e-3) * (cos(pi/rotf) - 1.0)
     &                    - (x(i,1)-dyn_org) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-1.025e-3) ) then ! middle tail
                disp(i,1) = -0.05 * (x(i,1)-dyn_org+2.125e-3) + xtsl
c     &                    + (x(i,1)-0.5e-3) * (cos(pi/rotf) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/rotf) - 1.0)
c     &                    + (x(i,1)-0.5e-3) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              else if ( x(i,1) .ge. (dyn_org+1.025e-3) ) then ! middle head
                disp(i,1) = -0.05 * (x(i,1)-dyn_org-2.125e-3) + xtsl
c     &                    + (x(i,1)-4.75e-3) * (cos(pi/rotf) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/rotf) - 1.0)
c     &                    + (x(i,1)-4.75e-3) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              else ! middle middle
                disp(i,1) = -0.05 * (x(i,1)-dyn_org) + xtsl
c     &                    + (x(i,1)-2.625e-3) * (cos(pi/rotf) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/rotf)
                disp(i,2) = -0.05 * (x(i,2)) ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/rotf) - 1.0)
c     &                    + (x(i,1)-2.625e-3) *  sin(pi/rotf)
                disp(i,3) = -0.05 * x(i,3)

              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            BC(i,1:3)   = disp(i,1:3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
      endif ! end if case 6
c
c.... end test case 6
c
c
c.... test case 7
c.... used to mov+rot+shrk ARO 9-grains case
c
      if ( casenumber .eq. 7 ) then
        dyn_org    = lstep * 0.04
        do i = 1,numnp
          if ( (ibits(iBC(i),14,3) .eq. 7) .and.
     &         (x(i,1) .lt. 23.9) .and. (x(i,1) .gt. -15.9) .and. 
     &         ((x(i,2)*x(i,2) + x(i,3)*x(i,3)) .lt. 63.9) ) then
            if( x(i,2) .ge. 2.0 ) then ! top
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! top tail
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm

              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! top head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! top middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2)-4.0)*(x(i,2)-4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
     &                    - (x(i,2)-4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)-4.0)/norm + 0.01
     &                    + (x(i,2)-4.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              endif ! end top--switch head middle tail

            else if( x(i,2) .le. -2.0 ) then ! bottom
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! bottom tail 
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
     &                    - (x(i,1)+8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! bottom head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
     &                    - (x(i,1)-8.0) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! bottom middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2)+4.0)*(x(i,2)+4.0)+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
     &                    - (x(i,1)) * (cos(pi/600) - 1.0)
     &                    + (x(i,2)+4.0) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2)+4.0)/norm - 0.01
     &                    - (x(i,2)+4.0) * (cos(pi/600) - 1.0)
     &                    - (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              endif ! end bottom--switch head middle tail

            else ! middle
              if ( x(i,1) .le. (dyn_org-4.0) ) then ! middle tail
                norm = SQRT((x(i,1)+8.0)*(x(i,1)+8.0)
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)+8.0)/norm + 0.02
c     &                    + (x(i,1)+8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)+8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else if ( x(i,1) .ge. (dyn_org+4.0) ) then ! middle head
                norm = SQRT((x(i,1)-8.0)*(x(i,1)-8.0)
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1)-8.0)/norm + 0.02
c     &                    + (x(i,1)-8.0) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)-8.0) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              else ! middle middle
                norm = SQRT((x(i,1))*(x(i,1))
     &                     +(x(i,2))*(x(i,2))+x(i,3)*x(i,3));
                disp(i,1) = -0.005 * (x(i,1))/norm + 0.02
c     &                    + (x(i,1)) * (cos(pi/600) - 1.0)
c     &                    - (x(i,2)) *  sin(pi/600)
                disp(i,2) = -0.005 * (x(i,2))/norm ! + 0.01
c     &                    + (x(i,2)) * (cos(pi/600) - 1.0)
c     &                    + (x(i,1)) *  sin(pi/600)
                disp(i,3) = -0.005 * x(i,3)/norm
 
              endif ! end middle--switch head middle tail

            endif ! end if switch top middle bottom
            BC(i,1:3)   = disp(i,1:3)
          endif ! end if inside cylinder
        enddo ! end loop numnp
      endif ! end if case 7
c
c.... end test case 7
c
      return
      end
c

