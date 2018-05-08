c----------------------------------------------------------------------
c
        module rigidBodyForce
          integer, allocatable :: rbIndex(:)
          real*8, allocatable  :: rbForce(:,:)
          real*8, allocatable  :: rbTorque(:,:)
        end module
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine malloc_rbForce ( )
c
        use rigidbody_m
        use rigidBodyForce
c
        allocate( rbForce(numrbs, 3)  )
        allocate( rbTorque(numrbs, 3) )
c
        return
        end
c
c----------------------------------------------------------------------
c
        subroutine init_rbForce ( )
c
        use rigidBodyForce
c
        rbForce  = 0.0
        rbTorque = 0.0
c
        return
        end
c
c----------------------------------------------------------------------
c
        subroutine commu_rbForce ( )
c
        use rigidBodyForce
        include "common.h"
c
        real*8, dimension(3) :: Forin, Forout
c
        do i = 1, numrbs
          if (iter .eq. nitr) then
            Forin  = (/ rbForce(i,1), rbForce(i,2), rbForce(i,3) /)
            if (numpe > 1) then
              call MPI_ALLREDUCE ( Forin(1),  Forout(1), 3,
     &                             MPI_DOUBLE_PRECISION,
     &                             MPI_SUM, MPI_COMM_WORLD, ierr)
              Force(i,1:3) = Forout(1:3)
            endif
          endif
        enddo
c
        return
        end
c
c----------------------------------------------------------------------
c
        subroutine release_rbForce ( )
c
        use rigidBodyForce
c
        deallocate( rbForce )
        deallocate( rbTorque )
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine local_rbIndex ( ien )
c
        use rigidBodyFlag
        use rigidBodyForce
c
        include "common.h"
c
        dimension    ien(npro,nshl)
        dimension    lnode(27)
        integer, dimension(npro, nshl) :: rbFlagsl
        integer   i,  j,  k,  counter
c
        allocate ( rbIndex(npro) )
        rbIndex = 0 ! 0 means it isn't a rigid body
c
c.... localize rigid body flag
c
        do i = 1, nshl
          rbFlagsl(:,i) = rbFlags(ien(:,i))
        enddo
c
c.... get local node index
c
        call getbnodes(lnode)
c
c.... if all boundary nodes share the same flag, set element flag
c
        do i = 1, npro
          counter = 0
          do j = 1, nenbl
            if ( (rbFlagsl(i, lnode(j)) .gt. 0) .and.
     &           (rbFlagsl(i, lnode(j)) .eq. rbFlagsl(i, lnode(1))) ) then
              counter = counter + 1
            endif
          enddo
          if (counter .eq. nenbl) then ! we find a mesh face on rigid body
            rbIndex(i) = rbFlagsl(i, lnode(1))
          endif
        enddo
c
        return
        end
c
c----------------------------------------------------------------------
c
        subroutine release_rbIndex ( )
c
        use rigidBodyForce
c
        deallocate( rbIndex )
c
        return
        end
c
c----------------------------------------------------------------------
c

c
c----------------------------------------------------------------------
c

c
c----------------------------------------------------------------------
c

c
c----------------------------------------------------------------------
c

