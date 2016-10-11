      module eqn_state_m
c
        use e3_def_m
        use number_def_m
        use matdat_def_m
c
        implicit none
c
        real*8, dimension(:), pointer :: rho,ei,p,T,h,cv,cp,alphaP,betaT,gamb,c
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
      subroutine getthm_solid_1
c
c
      end subroutine getthm_solid_1
c
      end module eqn_state_m
