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
          real*8,  dimension(npro,nsd), intent(in) :: nv
          real*8,  dimension(npro,nsd,nsd), intent(out) :: proj 
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
       end module e3if_dc_m
