      module eqn_state_m
c
        use e3_param_m
        use number_def_m
        use matdat_def_m
        use mmatpar_m
c
        implicit none
c
        real*8 :: rho_ref, p_ref, T_ref, alpha_P, beta_T, cv_liq
        real*8 :: rho_ref_s, p_ref_s, T_ref_s, alpha_P_s, cv_s
        real*8 :: bulkMod_s, shearMod_s 
c
      contains
c
      subroutine getthm6_ideal_gas
c
        mw    = mat_prop(mater,iprop_ideal_gas_mw, 1)
        gamma = mat_prop(mater,iprop_ideal_gas_gamma,1)
        Rgas  = 8314.0d0/mw
        gamma1 = gamma - one
c
        rho = pres / (Rgas*T)
        ei  = T * Rgas / gamma1
c
      end subroutine getthm6_ideal_gas
c
      subroutine getthm7_ideal_gas
c
        call getthm6_ideal_gas
c
        h   = T * Rgas / gamma1 * gamma
        cp  = Rgas*gamma / gamma1
        alphaP = one / T
        betaT  = one / pres
        if (associated(cv)) cv  = Rgas / gamma1
        if (associated(gamb)) gamb = gamma1
        if (associated(c)) c =  sqrt( gamma * Rgas * T )
c
      end subroutine getthm7_ideal_gas
c
      function rho_ideal_gas(p,R,T) result(rho)
        implicit none
        real*8 p,R,T,rho
        rho = p / (R*T)
      end function rho_ideal_gas
c
      subroutine getthm6_liquid_1
c
        rho_ref = mat_prop(mater,iprop_liquid_1_rho_ref,1)
        p_ref   = mat_prop(mater,iprop_liquid_1_p_ref,  1)
        T_ref   = mat_prop(mater,iprop_liquid_1_T_ref,  1)
        cv_liq  = mat_prop(mater,iprop_liquid_1_cv,     1)
        alpha_P = mat_prop(mater,iprop_liquid_1_alphaP, 1)
        beta_T  = mat_prop(mater,iprop_liquid_1_betaT,  1)
c
        rho = rho_ref * (one - alpha_P*(T-T_ref) + beta_T*(pres-P_ref))
        ei  = cv_liq*T
c
      end subroutine getthm6_liquid_1
c
      subroutine getthm7_liquid_1
c
        call getthm6_liquid_1
c
        h   = ei + pres/rho
        cp  = cv_liq
        alphaP = alpha_P
        betaT  = beta_T
        if (associated(cv)) cv  = cv_liq
        if (associated(gamb)) gamb = zero
c        c =  sqrt(one/(rho_ref*betaT))
        if (associated(c)) c =  sqrt(one/(rho*betaT))
c
      end subroutine getthm7_liquid_1
c
      subroutine getthm6_solid_1
c
        use e3_solid_m
        implicit none
c
        rho_ref_s = mat_prop(mater,iprop_solid_1_rho_ref,1)
        p_ref_s   = mat_prop(mater,iprop_solid_1_p_ref,  1)
        T_ref_s   = mat_prop(mater,iprop_solid_1_T_ref,  1)
        cv_s     = mat_prop(mater,iprop_solid_1_cv,     1)
        alpha_P_s = mat_prop(mater,iprop_solid_1_alphaP, 1)
        bulkMod_s = mat_prop(mater,iprop_solid_1_bulkMod, 1)
        shearMod_s = mat_prop(mater,iprop_solid_1_shearMod, 1)
c        beta_T  = mat_prop(mater,iprop_solid_1_betaT,  1)
c
        alphaP = alpha_P_s
        betaT  = one /(bulkMod_s * Ja_def) ! double check here
        rho = rho_ref_s * (one - alphaP*(T - T_ref_s) 
     &                    + betaT*(pres - p_ref_s))
c        ei  = ( cv_s - pres * alpha_P_s/rho)* T 
c    &        + (betaT * pres - alphaP * T)/rho * pres
        ei  = cv_s * T
c
      end subroutine getthm6_solid_1
c
      subroutine getthm7_solid_1
c
        use e3_solid_m
        implicit none
c
        call getthm6_solid_1
c
        h   = ei + pres/rho
        cv  = cv_s
        cp  = cv_s
        bulkMod = bulkMod_s
        shearMod = shearMod_s
c        c =  sqrt(one/(rho_ref*betaT))
C         c =  sqrt(one/(rho*betaT))
C         gamb = zero
c
c
      end subroutine getthm7_solid_1
c
      subroutine getthmif
c
        use e3if_param_m
c
        rho  => rho0
        ei   => ei0
        pres => pres0
        T    => T0
        h    => h0
        cp   => cp0
        alphaP => alfaP0
        betaT => betaT0
        mater = mater0
c
        call getthmif0_ptr
c
        rho  => rho1
        ei   => ei1
        pres => pres1
        T    => T1
        h    => h1
        cp   => cp1
        alphaP => alfaP1
        betaT => betaT1
        mater = mater1
c
        call getthmif1_ptr
c
      end subroutine getthmif
c
      end module eqn_state_m
