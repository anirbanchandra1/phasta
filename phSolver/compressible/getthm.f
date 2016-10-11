      subroutine getthm (rho_,ei_,p_,T_,npro_,mater_
     &,                  h_,  cv_,cp_,alphaP_,betaT_,gamb_,c_)
        use eqn_state_m
c
        implicit none
c
        integer, intent(in) :: npro_, mater_
        real*8, dimension(npro), target, intent(in) :: p_,T_
        real*8, dimension(npro), target, intent(inout) :: rho_,ei_,h_,cv_,cp_,alphaP_,betaT_,gamb_,c_
c
        npro = npro_
        mater = mater_
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
        gamb  => gamb_
        c => c_
c 
        select case (mat_eos(mater,1))
        case (ieos_ideal_gas,ieos_ideal_gas_2)
          call getthm_ideal_gas
        case (ieos_liquid_1)
          call getthm_liquid_1
        case default
          call error ('getthm  ', 'wrong material', mater)
        end select
c
      end subroutine getthm
