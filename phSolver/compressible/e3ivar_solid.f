      subroutine e3ivar_solid(rho_,     ei_,      p_,       T_,      h_,
     &                        cv_,      cp_,      alphaP_,   betaT_,
     &                        bulkMod_, shearMod_,Ja_def_,   d_,
     &                        det_d_,   det_baf_, g1yi_,    g2yi_,    g3yi_,
     &                        npro_,    nsd_,     almBi_,    alfBi_,
     &                        gamBi_,    intp_,   Delt_)
c
        use solid_func_m
c
        implicit none
c
        real*8, dimension(npro_), intent(out) :: bulkMod_, shearMod_
        real*8, dimension(npro_), target, intent(inout) :: rho_, ei_,p_,
     &                                                     T_, h_, cv_,
     &                                                     cp_, alphaP_,
     &                                                     betaT_
        real*8, dimension(npro_,nsd_), target, intent(in) :: g1yi_,g2yi_,g3yi_
        real*8, dimension(npro_,6), target :: d_
        real*8, dimension(npro_), target :: Ja_def_
        real*8, dimension(npro_), target :: det_d_
        real*8, dimension(npro_), target :: det_baf_
        real*8, intent(in) :: Delt_
        integer, intent(in) :: npro_, nsd_
        integer, intent(in) :: intp_
        real*8, intent(in) :: almBi_, alfBi_, gamBi_
c
        npro = npro_
        nsd = nsd_
        Delt_s = Delt_
        intp_s = intp_
        almBi_s = almBi_
        alfBi_s = alfBi_
        gamBi_s = gamBi_  
c
        rho => rho_
        ei  => ei_
        p   => p_
        T   => T_
        h   => h_
        cv  => cv_
        cp  => cp_
        alphaP => alphaP_
        betaT  => betaT_
c        gamb  => gamb_
c        c => c_
c        bulkMod =>bulkMod_
c        shearMod =>shearMod_
c
        d => d_
        Ja_def => Ja_def_
        det_d => det_d_
        det_baf => det_baf_
        dudx => g1yi_(:,2:4)
        dudy => g2yi_(:,2:4)
        dudz => g3yi_(:,2:4)
c
        call calc_solid
        call getthm_solid_1(bulkMod_, shearMod_,ja_def_)
c
      end subroutine e3ivar_solid
