      module core_snap
      use iso_c_binding
c
      interface
c
        subroutine core_get_pos_on_surf ( dx, dy, dz, numnp, px, py, pz )
     &    bind(C, NAME='core_get_pos_on_surf')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(in),dimension(:) ::
     &                               dx(numnp), dy(numnp), dz(numnp)
          real(c_double),intent(out),dimension(:) ::
     &                               px(numnp), py(numnp), pz(numnp)
        end subroutine
c
        subroutine core_is_in_closure(e_dim, e_tag, t_dim, t_tag, answer)
     &    bind(C, NAME='core_is_in_closure')
        use iso_c_binding
          integer(c_int), intent(in) :: e_dim, e_tag, t_dim, t_tag
          integer(c_int), intent(out) :: answer
        end subroutine
c
      end interface
c
      end module
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      module core_mesh_quality
      use iso_c_binding
c
      interface
c
        subroutine core_measure_mesh(x1, x2, x3, numnp, minvq, minfq)
     &    bind(C, NAME='core_measure_mesh')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(in),dimension(:) ::
     &                               x1(numnp), x2(numnp), x3(numnp)
          real(c_double),intent(out) :: minvq, minfq
        end subroutine
c
      end interface
c
      end module
