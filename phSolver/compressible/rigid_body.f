c----------------------------------------------------------------------
c
        module rigidBodyForce
          integer, allocatable :: rbIndex(:)
          real*8, allocatable  :: rbForce(:,:)
          real*8, allocatable  :: rbTorque(:,:)
c
          real*8, allocatable  :: rbDisp(:,:)
          real*8, allocatable  :: rbTotalDisp(:,:)
          real*8, allocatable  :: rbForceOld(:,:)
          real*8, allocatable  :: rbTorqueOld(:,:)
          real*8, allocatable  :: rbVelOld(:,:)
          real*8, allocatable  :: rbAccOld(:,:)
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
        allocate( rbDisp(numrbs, 3) )
        allocate( rbTotalDisp(numrbs, 3) )
        allocate( rbForceOld(numrbs, 3)  )
        allocate( rbTorqueOld(numrbs, 3)  )
        allocate( rbVelOld(numrbs, 3)  )
        allocate( rbAccOld(numrbs, 3)  )
c
        rbForce  = 0.0
        rbTorque = 0.0
c
        rbDisp = 0.0
        rbTotalDisp = 0.0
        rbForceOld = 0.0
        rbTorqueOld = 0.0
        rbVelOld = 0.0
        rbAccOld = 0.0
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
        use rigidBodyFlag
        use rigidBodyForce
c
        include "common.h"
        include "mpif.h"
c
        real*8, dimension(3) :: Forin, Forout
c
        do i = 1, numrbs
c.... XXX need to add rotation force
          if (rbsMM(i) .ne. 1)
     &      call error('rigidBodyBCElas','not support mode',rbsMM(i))
          if (iter .eq. nitr) then
            Forin  = (/ rbForce(i,1), rbForce(i,2), rbForce(i,3) /)
            if (numpe > 1) then
              call MPI_ALLREDUCE ( Forin,  Forout, 3,
     &                             MPI_DOUBLE_PRECISION,
     &                             MPI_SUM, MPI_COMM_WORLD, ierr)
              rbForce(i,1:3) = Forout(1:3)
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
        deallocate( rbDisp )
        deallocate( rbTotalDisp )
        deallocate( rbForceOld )
        deallocate( rbTorqueOld )
        deallocate( rbVelOld )
        deallocate( rbAccOld )
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
        subroutine calc_rbMotion ( )
c
        use rigidBodyFlag
        use rigidBodyForce
c
        include "common.h"
c
        real*8                           t_mag
        real*8, dimension(nsd)        :: t_dir,  tmpAcc
        real*8, dimension(numrbs,nsd) :: t_acc
c.... generalized-alpha parameters
        real*8  am, af, gm, bt
c
        am = 0.5 * (3.0-rhoinf_rb(1)) / (1.0+rhoinf_rb(1))
        af = 1.0 / (1.0+rhoinf_rb(1))
        gm = 0.5 + am - af
        bt = 0.25 * (1.0+am-af) * (1.0+am-af)
c
c.... loop over rigid body
        do j = 1,numrbs
          if (rbsMM(j) .ne. 1)
     &      call error('rigidBodyBCElas','not support mode',rbsMM(j))
c.... debugging {
c.... XXX need to implement rotation motion
c          call core_get_centroid(rbMTs(j), cent(j,:))
c.... debugging }
c
c.... apply constraint on force
c
          t_dir(:) = rb_prop(j,2:4) ! translation direction
          t_mag = sqrt( t_dir(1)*t_dir(1)
     &                + t_dir(2)*t_dir(2)
     &                + t_dir(3)*t_dir(3) )
          if (t_mag .lt. 1e-12) then ! no constraint
            t_dir(1:3) = 1.0
          else
            t_dir(1:3) = t_dir(1:3) / t_mag
          endif
c
          rbForce(j,1:3) = rbForce(j,1:3) * t_dir(1:3)
c.... debugging {
c          rbForce(j,1) = 3.0 + 2.0 * Delt(1) * (lstep + 1.0)
c.... debugging }
c
c.... set initial acceleration
c
          if (istep .eq. 0) then
            rbAccOld(j,1:3) = rbForce(j,1:3) / rb_prop(j,1)
            rbForceOld(j,1:3) = rbForce(j,1:3)
          endif
c
c.... prepare variable
c
          tmpAcc(1:3) = af/am/rb_prop(j,1) * rbForce(j,1:3)
     &                + (1-af)/am/rb_prop(j,1) * rbForceOld(j,1:3)
c
c.... get displacement
c
          rbDisp(j,1:3) = Delt(1) * rbVelOld(j,1:3)
     &                  + Delt(1) * Delt(1) * (
     &                    (0.5 + bt/am) * rbAccOld(j,1:3)
     &                  - bt * tmpAcc(1:3)    )
c
c.... get total displacement
c
          rbTotalDisp(j,1:3) = rbTotalDisp(j,1:3) + rbDisp(j,1:3)
c
c.... get velocity
c
          rbVelOld(j,1:3) = rbVelOld(j,1:3)
     &                    + Delt(1) * (1.0 - gm/am) * rbAccOld(j,1:3)
     &                    + Delt(1) * gm * tmpAcc(1:3)
c
c.... get acceleration
c
          rbAccOld(j,1:3) = tmpAcc(1:3) - (1.0-am)/am * rbAccOld(j,1:3)
c
c.... debugging {
c          if (myrank .eq. master)
c     &      write(*,*) "rbForce", rbForce(j,1)
c     &                ,"disp:", rbDisp(j,1)
c     &                ,"Velo:", rbVelOld(j,1)
c.... debugging }
c
c.... update force of the previous time step
c
          rbForceOld(j,1:3) = rbForce(j,1:3)
c
        enddo
c
        return
        end
c
c----------------------------------------------------------------------
c
        subroutine synchronize_rbForce ( )
c
        use rigidBodyFlag
        use rigidBodyForce
        use core_rigid_body
c
        include "common.h"
c
        real*8, dimension(numrbs) :: ax, ay, az, px, py, pz, ag, sc
c
c.... XXX need to add rotation motion here
        do j = 1,numrbs
          if (rbsMM(j) .ne. 1)
     &      call error('rigidBodyBCElas','not support mode',rbsMM(j))
        enddo
c
        ax = 1.0
        ay = 0.0
        az = 0.0
        px = 0.0
        py = 0.0
        pz = 0.0
        ag = 0.0
        sc = 0.0
c
        call core_update_rbms(rbTotalDisp(:,1),
     &                        rbTotalDisp(:,2),
     &                        rbTotalDisp(:,3),
     &                        ax, ay, az,
     &                        px, py, pz,
     &                        ag, sc,
     &                        rbsTags(:), numrbs)
c
        return
        end
c
c
c----------------------------------------------------------------------
c

c
c----------------------------------------------------------------------
c

