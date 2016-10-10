      subroutine e3ivar_solid(g1yi_,g2yi_,g3yi_,rho_,npro_,nsd_)
c
        use solid_func_m
c
        implicit none
c
        real*8, dimension(npro_,nsd_), target, intent(in) :: g1yi_,g2yi_,g3yi_
        real*8, dimension(npro_), target, intent(inout) :: rho_
        integer, intent(in) :: npro_,nsd_
c
        npro = npro_
        nsd = nsd_
c
        dudx => g1yi_(:,2:4)
        dudy => g2yi_(:,2:4)
        dudz => g3yi_(:,2:4)
c
        call calc_solid
        call getthm_solid_1(rho_,npro_)
c
      end subroutine e3ivar_solid
