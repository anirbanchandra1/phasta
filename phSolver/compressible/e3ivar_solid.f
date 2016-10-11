      subroutine e3ivar_solid(rho_, bulkMod_, shearMod_, Ja_def_, g1yi_,g2yi_,g3yi_,npro_,nsd_)
c
        use solid_func_m
c
        implicit none
c
        real*8, dimension(npro_), intent(out) :: bulkMod_, shearMod_
        real*8, dimension(npro_), target, intent(inout) :: rho_
        real*8, dimension(npro_,nsd_), target, intent(in) :: g1yi_,g2yi_,g3yi_
        real*8, dimension(npro_,nsd_), intent(in) :: Ja_def_
        integer, intent(in) :: npro_,nsd_
c
        npro = npro_
        nsd = nsd_
c
        rho => rho_
c        ei  => ei_
c        p   => p_
c        T   => T_
c        h   => h_
c        cv  => cv_
c        cp  => cp_
c        alphaP => alphaP_
c        betaT  => betaT_
c        gamb  => gamb_
c        c => c_
c
        dudx => g1yi_(:,2:4)
        dudy => g2yi_(:,2:4)
        dudz => g3yi_(:,2:4)
c
        call calc_solid
        call getthm_solid_1(bulkMod_, shearMod_, Ja_def_)
c
      end subroutine e3ivar_solid
