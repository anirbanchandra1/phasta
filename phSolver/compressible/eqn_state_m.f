      module eqn_state_m
c
        use e3_def_m
        use number_def_m
        use matdat_def_m
c
        implicit none
c
        real*8, dimension(:), pointer :: rho, ei, p, T, h, cv, cp,
     &                                   alphaP, betaT, gamb, c
c
      contains
c
      subroutine getthm_ideal_gas
c
        real*8 :: Rgas,gamma,gamma1,mw
c
        mw    = mat_prop(mater,iprop_ideal_gas_mw, 1)
        gamma = mat_prop(mater,iprop_ideal_gas_gamma,1)
        Rgas  = 8314.0d0/mw
        gamma1 = gamma - one
c
        rho = P / (Rgas*T)
        ei  = T * Rgas / gamma1
        h   = T * Rgas / gamma1 * gamma
        cv  = Rgas / gamma1
        cp  = Rgas*gamma / gamma1
        alphaP = one / T
        betaT  = one / P
        gamb = gamma1
        c =  sqrt( gamma * Rgas * T )
c
      end subroutine getthm_ideal_gas
c
      function rho_ideal_gas(p,R,T) result(rho)
        implicit none
        real*8 p,R,T,rho
        rho = p / (R*T)
      end function rho_ideal_gas
c
      subroutine getthm_liquid_1
c
        real*8 :: rho_ref, p_ref, T_ref, alpha_P, beta_T, cv_
c
        rho_ref = mat_prop(mater,iprop_liquid_1_rho_ref,1)
        p_ref   = mat_prop(mater,iprop_liquid_1_p_ref,  1)
        T_ref   = mat_prop(mater,iprop_liquid_1_T_ref,  1)
        cv_     = mat_prop(mater,iprop_liquid_1_cv,     1)
        alpha_P = mat_prop(mater,iprop_liquid_1_alphaP, 1)
        beta_T  = mat_prop(mater,iprop_liquid_1_betaT,  1)
c
        rho = rho_ref * (one - alpha_P*(T-T_ref) + beta_T*(P-P_ref))
        ei  = cv_*T
        h   = ei + P/rho
        cv  = cv_
        cp  = cv_
        alphaP = alpha_P
        betaT  = beta_T
c        c =  sqrt(one/(rho_ref*betaT))
        c =  sqrt(one/(rho*betaT))
        gamb = zero
c
      end subroutine getthm_liquid_1
c
c
      subroutine getthm_solid_1(bulkMod, shearMod,Ja_def)
c
        real*8, dimension(npro), intent(out) :: bulkMod, shearMod
        real*8, dimension(npro), intent(in)  :: Ja_def
c
        real*8 :: rho_ref_s, p_ref_s, T_ref_s, alpha_P_s, cv_s
        real*8 :: bulkMod_s, shearMod_s 
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
     &                    + betaT*(p - p_ref_s))
c        ei  = ( cv_s - P * alpha_P_s/rho)* T 
c    &        + (betaT * P - alphaP * T)/rho * p
        ei  = cv_s * T
        h   = ei + p/rho
        cv  = cv_s
        cp  = cv_s
        bulkMod = bulkMod_s
        shearMod = shearMod_s
c        c =  sqrt(one/(rho_ref*betaT))
C         c =  sqrt(one/(rho*betaT))
C         gamb = zero
c
c
      end subroutine getthm_solid_1
c
      end module eqn_state_m
