      module e3if_lhs_m
c
        use e3if_defs_m
c
        implicit none
c
      contains
c
      subroutine set_lhs_matrices
c
        integer :: isd,jsd,n,p,q,r
c
        do q = 1,nflow
          do p = 1,nflow
            do isd = 1,nsd
c
              AiNa0  (:,isd,p,q) = zero
              KijNaj0(:,isd,p,q) = zero
c
              do n = 1,nshl0
                AiNa0  (:,isd,p,q) = AiNa0(:,isd,p,q) + Ai0(:,isd,p,q)*shp0(:,n)
                do jsd = 1,nsd
                  KijNaj0(:,isd,p,q) = KijNaj0(:,isd,p,q) + Kij0(:,isd,jsd,p,q)*shg0(:,n,jsd)
                enddo
              enddo
c
              AiNa1  (:,isd,p,q) = zero
              KijNaj1(:,isd,p,q) = zero
c
              do n = 1,nshl1
                AiNa1  (:,isd,p,q) = AiNa1(:,isd,p,q) + Ai1(:,isd,p,q)*shp1(:,n)
                do jsd = 1,nsd
                  KijNaj1(:,isd,p,q) = KijNaj1(:,isd,p,q) + Kij1(:,isd,jsd,p,q)*shg1(:,n,jsd)
                enddo
              enddo
c
            enddo
          enddo
        enddo
c
        do q = 1,nflow
          do p = 1,nflow
            do isd = 1,nsd
c
              KijNajC0(:,isd,p,q) = zero
              KijNajC1(:,isd,p,q) = zero
c
              do r = 1,nflow
                KijNajC0(:,isd,p,q) = KijNajC0(:,isd,p,q) + KijNaj0(:,isd,p,r)*cmtrx(:,r,q)
                KijNajC1(:,isd,p,q) = KijNajC1(:,isd,p,q) + KijNaj1(:,isd,p,r)*cmtrx(:,r,q)
              enddo
c
            enddo
          enddo
        enddo
c
        return
c
      end subroutine set_lhs_matrices
c
      subroutine calc_egmass( egmass_self_self,    egmass_self_nbr,
     &                         AiNa_self, AiNa_nbr, KijNaj_self, KijNaj_nbr,
     &                         shp, n, WdetJ,
     &                         nshl_self, nshl_nbr)
c        
        real*8, dimension(:,:,:), intent(inout) :: egmass_self_self, egmass_self_nbr
        real*8, dimension(:,:,:,:), intent(in)  :: AiNa_self, AiNa_nbr, KijNaj_self, KijNaj_nbr
        real*8, dimension(:,:), intent(in) :: shp, n
        real*8, dimension(:), intent(in) :: WdetJ
        integer, intent(in) :: nshl_self,nshl_nbr
c
        integer :: i,j,p,q,inode0,inode1,isd
        integer :: i0,j0,il,jl,iflow,jflow
c
        do i = 1,nshl_self
          i0 = nflow*(i-1)
c
          do j = 1,nshl_self
            j0 = nflow*(j-1)
c
            do iflow = 1,nflow
              do jflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                do isd = 1,nsd
                  egmass_self_self(:,il,jl) = egmass_self_self(:,il,jl) - (
     &              pt50 * shp(:,i) * ( AiNa_self(:,isd,iflow,jflow) 
     &                                - KijNaj_self(:,isd,iflow,jflow) ) * n(:,isd) 
     &              ) * WdetJ
                enddo
c
              enddo
            enddo
c
          enddo
c
          do j = 1,nshl_nbr
            j0 = nflow*(j-1)
c
            do iflow = 1,nflow
              do jflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                do isd = 1,nsd
                  egmass_self_nbr(:,il,jl) = egmass_self_nbr(:,il,jl) - (
     &              pt50 * shp(:,i) * ( AiNa_nbr(:,isd,iflow,jflow) 
     &                                - KijNaj_nbr(:,isd,iflow,jflow) ) * n(:,isd) 
     &              ) * WdetJ
                enddo
c
              enddo
            enddo
c
          enddo
c
        enddo
c
      end subroutine calc_egmass
c
      subroutine calc_egmass_(egmass,AiNa1,KijNaj0,KijNaj1,KijNajC0,shp0,shp1,n0,n1,WdetJ,nshl0,nshl1)
c
        real*8, dimension(:,:,:), intent(inout) :: egmass
        real*8, dimension(:,:,:,:), intent(in) :: AiNa1,KijNaj0,KijNaj1,KijNajC0
        real*8, dimension(:,:), intent(in) :: shp0,shp1,n0,n1
        real*8, dimension(:), intent(in) :: WdetJ
        integer, intent(in) :: nshl0,nshl1
c
        integer :: i,j,p,q,inode0,inode1,isd
        integer :: i0,j0,il,jl,iflow,jflow
c
c... loop through the columns
c
        do j = 1, nshl1
          j0 = nflow*(j - 1)
c
          do i = 1,nshl0
            i0 = nflow*(i - 1)
c
            do jflow = 1,nflow
              do iflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                do isd = 1,nsd
c
                  egmass(:,il,jl) = egmass(:,il,jl) + (
     &                pt50 * shp0(:,i) * ( 
     &                AiNa1  (:,isd,iflow,jflow) - KijNaj1(:,isd,iflow,jflow)) * n0(:,isd)
     &             + pt50 * s * KijNajC0(:,isd,iflow,jflow)*n1(:,isd)
     &             + e*mu/h * ctc(:,iflow,jflow)*shp0(:,i)*n0(:,isd)*shp1(:,j)*n1(:,isd)
     &            ) * WdetJ
c
                enddo
              enddo
            enddo
          enddo
        enddo
c              
      return
c
      end subroutine calc_egmass_
c
      subroutine calc_egmass_2(egmass,AiNa1,KijNaj0,KijNaj1,KijNajC0,shp0,shp1,n0,n1,WdetJ,nshl0,nshl1)
c
        real*8, dimension(:,:,:), intent(inout) :: egmass
        real*8, dimension(:,:,:,:), intent(in) :: AiNa1,KijNaj0,KijNaj1,KijNajC0
        real*8, dimension(:,:), intent(in) :: shp0,shp1,n0,n1
        real*8, dimension(:), intent(in) :: WdetJ
        integer, intent(in) :: nshl0,nshl1
c
        integer :: i,j,p,q,inode0,inode1,isd
        integer :: i0,j0,il,jl,iflow,jflow
c
c... loop through the columns
c
        do j = 1, nshl1
          j0 = nflow*(j - 1)
c
          do i = 1,nshl0
            i0 = nflow*(i - 1)
c
            do jflow = 1,nflow
              do iflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                do isd = 1,nsd
c
                  egmass(:,il,jl) = egmass(:,il,jl) + (
     &              - pt50 * shp0(:,i) * ( 
c     &                AiNa1  (:,isd,iflow,jflow)
     &              + KijNaj1(:,isd,iflow,jflow)
     &              ) * n0(:,isd)
c     &             + pt50 * s * KijNajC0(:,isd,iflow,jflow)*n1(:,isd)
c     &             + e*mu/h * ctc(:,iflow,jflow)*shp0(:,i)*n0(:,isd)*shp1(:,j)*n1(:,isd)
     &            ) * WdetJ
c
                enddo
              enddo
            enddo
          enddo
        enddo
c              
      return
c
      end subroutine calc_egmass_2
c
      subroutine calc_egmass_old(egmass,AiNa1,KijNaj0,KijNaj1,KijNajC0,shp0,shp1,n0,n1,WdetJ,nshl0,nshl1)
c
        real*8, dimension(:,:,:), intent(inout) :: egmass
        real*8, dimension(:,:,:,:), intent(in) :: AiNa1,KijNaj0,KijNaj1,KijNajC0
        real*8, dimension(:,:), intent(in) :: shp0,shp1,n0,n1
        real*8, dimension(:), intent(in) :: WdetJ
        integer, intent(in) :: nshl0,nshl1
c
        integer :: i,j,p,q,inode0,inode1,isd
        integer :: i0,j0,il,jl,iflow,jflow
c
        do inode1 = 1, nshl1
          q = nflow * (inode1 - 1) 
c
          do inode0 = 1, nshl0
            p = nflow * (inode0 - 1)
c
            do isd = 1,nsd
              egmass(:,p,q) = egmass(:,p,q) + (
     &              pt50 * shp0(:,inode0) * (AiNa1(:,isd,p,q)+KijNaj1(:,isd,p,q))*n0(:,isd)
     &            + pt50 * s * KijNajC0(:,isd,p,q)*n1(:,isd)
     &            + e*mu/h * ctc(:,p,q)*shp0(:,inode0)*n0(:,isd)*shp1(:,inode1)*n1(:,isd)
     &            )*WdetJ
            enddo
c
          enddo
c
        enddo
c
      end subroutine calc_egmass_old
c
      end module e3if_lhs_m
