c
c  rho                   : density
c  p                     : pressure
c  T                     : temperature
c  ei                    : internal energy
c  h                     : enthalpy
c  alphaP                : expansivity
c  betaT                 : isothermal compressibility
c  cp                    : specific heat at constant pressure
c  rk                    : kinetic energy
c
      module e3_param_m
c
        use propar_m
c
        implicit none
c
        abstract interface
          subroutine getthm2
          end subroutine getthm2
        end interface
c
        integer :: nqpt
        integer :: mater0, mater1, mater
c 
        real*8, dimension(:), pointer :: rho, ei, p, pres, T, h, cv, cp,
     &                             alphaP, betaT, gamb, c, rk
c
        procedure(e3_malloc), pointer :: e3_malloc_ptr
        procedure(e3_mfree), pointer :: e3_mfree_ptr
        procedure(getthm2), pointer :: getthm6_ptr, getthm7_ptr
c
        contains
c
        subroutine e3_malloc
          allocate(rho(npro))
          allocate(ei(npro))
          allocate(pres(npro))
          allocate(p(npro))
          allocate(T(npro))
          allocate(h(npro))
          allocate(cv(npro))
          allocate(cp(npro))
          allocate(alphap(npro))
          allocate(betat(npro))
          allocate(gamb(npro))
          allocate(c(npro))
          allocate(rk(npro))
        end subroutine e3_malloc
c
        subroutine e3_mfree
          deallocate(rho)
          deallocate(ei)
          deallocate(pres)
          deallocate(p)
          deallocate(T)
          deallocate(h)
          deallocate(cv)
          deallocate(cp)
          deallocate(alphap)
          deallocate(betat)
          deallocate(gamb)
          deallocate(c)
          deallocate(rk)
        end subroutine e3_mfree
c
      end module e3_param_m
