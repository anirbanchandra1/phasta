      module e3if_dc_m
c
c------------------------------------------------------------------------------
c  calculating the discontinuity capturing (DC) operator for the interface
c  macro-elements
c------------------------------------------------------------------------------
        implicit none
c
c
      contains
        subroutine calc_projector(proj, nv)
c..............................................................................
c  calculating the projector which project a quanlity into its tangential
c  components as :
c  \bm{P} =  \bm{I} - \bm{n} \otimes \bm{n}
c..............................................................................        
          use number_def_m
          use propar_m, only: npro
          use global_const_m, only: nsd
          implicit none
c          
          real*8,  dimension(npro,nsd), intent(in) :: nv ! normal vector
          real*8,  dimension(npro,nsd,nsd), intent(out) :: proj ! projector
c
          proj = zero
c
          proj(:,1,1) = one - nv(:,1)*nv(:,1)
          proj(:,1,2) =     - nv(:,1)*nv(:,2)
          proj(:,1,3) =     - nv(:,1)*nv(:,3)
c
          proj(:,2,1) =     - nv(:,2)*nv(:,1)
          proj(:,2,2) = one - nv(:,2)*nv(:,2)
          proj(:,2,3) =     - nv(:,2)*nv(:,3)
c
          proj(:,3,1) =     - nv(:,3)*nv(:,1)
          proj(:,3,2) =     - nv(:,3)*nv(:,2)
          proj(:,3,3) = one - nv(:,3)*nv(:,3)
  
c
        end subroutine calc_projector
c
        subroutine calc_ch(ch0, ch1, f_jump)
c..............................................................................
c  calculating the c^h of the DC operator for the interface, which is analogous 
c  to the nu^h for the volumn elements. 
c  Unit of c^h ~ [Length]/ [Time]
c..............................................................................        
          use number_def_m
          use propar_m, only: npro
          use conpar_m, only: nflow
          implicit none
c
          real*8, dimension(npro),intent(out) :: ch0, ch1 
          real*8, dimension(npro,nflow),intent(in) :: f_jump ! flux jump
c
          real*8, dimension(npro,nflow) :: u_ref_0, u_ref_1 ! reference conservative
                                                            ! variables
          real*8, dimension(npro) :: factor ! non-dimensional parameter
          real*8, dimension(npro,nflow) :: temp0, temp1 ! local temporary array
          integer :: iel                                                  
c
          factor = one
c
          u_ref_0(:,1) = 1.0d0 ! hacking, rho_ref, gas phase
          u_ref_0(:,2) = 1.0d0 * 1.0d2 ! hacking, rho_ref*v_ref
          u_ref_0(:,3) = 1.0d0 * 1.0d2 ! hacking, rho_ref*v_ref
          u_ref_0(:,4) = 1.0d0 * 1.0d2 ! hacking, rho_ref*v_ref
          u_ref_0(:,5) = 1.0d0 * (1.0d2)**two ! hacking, rho_ref*v_ref^2
c
          u_ref_1(:,1) = 1.0d2 ! hacking, rho_ref, gas phase
          u_ref_1(:,2) = 1.0d2 * 1.0d0 ! hacking, rho_ref*v_ref
          u_ref_1(:,3) = 1.0d2 * 1.0d0 ! hacking, rho_ref*v_ref
          u_ref_1(:,4) = 1.0d2 * 1.0d0 ! hacking, rho_ref*v_ref
          u_ref_1(:,5) = 1.0d2 * (1.0d0)**two ! hacking, rho_ref*v_ref^2
c... diag(U1_ref ... U5_ref) * flux_jump          
          temp0(:,1) = ( one/u_ref_0 (:,1) )*f_jump(:,1)
          temp0(:,2) = ( one/u_ref_0 (:,2) )*f_jump(:,2)
          temp0(:,3) = ( one/u_ref_0 (:,3) )*f_jump(:,3)
          temp0(:,4) = ( one/u_ref_0 (:,4) )*f_jump(:,4)
          temp0(:,5) = ( one/u_ref_0 (:,5) )*f_jump(:,5)
c
          temp1(:,1) = ( one/u_ref_1 (:,1) )*f_jump(:,1)
          temp1(:,2) = ( one/u_ref_1 (:,2) )*f_jump(:,2)
          temp1(:,3) = ( one/u_ref_1 (:,3) )*f_jump(:,3)
          temp1(:,4) = ( one/u_ref_1 (:,4) )*f_jump(:,4)
          temp1(:,5) = ( one/u_ref_1 (:,5) )*f_jump(:,5)
c... c^h for each phase
          do iel = 1,npro
            ch0(iel) = factor(iel) 
     &               * sqrt(dot_product(temp0(iel,:),temp0(iel,:)))
            ch1(iel) = factor(iel) 
     &               * sqrt(dot_product(temp1(iel,:),temp1(iel,:)))
          enddo
c          
        end subroutine calc_ch
c
c
      end module e3if_dc_m
