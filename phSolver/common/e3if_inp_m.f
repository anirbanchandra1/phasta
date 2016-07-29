      module e3if_inp_m
c
        implicit none
c
        integer :: vi_ramping, vi_model
        real*8 :: ramp_time
        real*8 :: vi_mag, dgif_alpha, dgif_beta, dgif_s, dgif_e, dgif_h
c
        common /e3if_dat/ vi_ramping, vi_model, ramp_time, vi_mag, dgif_alpha, dgif_beta, dgif_s, dgif_e, dgif_h
c
        integer, parameter :: no_ramp = 1, 
     &                        linear_ramp = 2
c
        integer, parameter :: const_vi = 1,
     &                        vi_model1 = 2
c
      end module e3if_inp_m
