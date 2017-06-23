      module e3if_lhs_m
c
c----------------------------------------
c    detail subroutines to computer the
c    LHS for interface elements
c----------------------------------------
        use e3if_param_m
        use e3if_func_m
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
c
            A0Na_0(:,p,q) = zero
            A0Na_1(:,p,q) = zero
c
            do n = 1,nshl0
              A0Na_0(:,p,q) = A0Na_0(:,p,q) + A0_0(:,p,q)*shp0(:,n)
            enddo
c
            do n = 1,nshl1
              A0Na_1(:,p,q) = A0Na_1(:,p,q) + A0_1(:,p,q)*shp1(:,n)
            enddo
c
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
      subroutine calc_egmass_ (
     &                        egmass
     &,                       shp0, shp1, shg0, shg1
     &,                       Ai1
     &,                       Kij0, Kij1
     &,                       ni0, ni1, WdetJ0
     &,                       nshl0, nshl1
     &                        )
        real*8, dimension(:,:,:),   intent(inout) :: egmass
        real*8, dimension(:,:),     intent(in)    :: shp0, shp1
        real*8, dimension(:,:,:),   intent(in)    :: shg0, shg1
        real*8, dimension(:,:,:,:), intent(in)    :: Ai1
        real*8, dimension(:,:,:,:,:), intent(in) :: Kij0, Kij1
        real*8, dimension(:,:),     intent(in)    :: ni0, ni1
        real*8, dimension(:),       intent(in)    :: WdetJ0
        integer, intent(in)  :: nshl0,nshl1
c
        integer :: i,j,i0,j0,il,jl,iflow,jflow,kflow,isd,jsd
        real*8 :: this_mu(npro,nflow,nflow)
        real*8, dimension(npro) :: Ai1Na1ni0,Kij1Naj1ni0,Kij0Nbj0CNa1ni1,CNb0ni0CNa1ni1
c
c---------> ATTENTION <--------
c  This still needs to be confirmed...
c
      this_mu = zero
      this_mu(:,2,2) = mu(:,2) ! mu
      this_mu(:,3,3) = mu(:,3) ! "
      this_mu(:,4,4) = mu(:,4) ! "
      this_mu(:,5,5) = mu(:,5) ! kappa
c
c---------> END ATTANETION <-----
c
        do i = 1,nshl0
          i0 = nflow*(i-1)
c
          do j = 1,nshl0
            j0 = nflow*(j-1)
c
            do iflow = 1,nflow
              do jflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                Ai1Na1ni0       = zero
                Kij1Naj1ni0     = zero
                Kij0Nbj0CNa1ni1 = zero
                CNb0ni0CNa1ni1  = zero
c
                do isd = 1,nsd
                  Ai1Na1ni0 = Ai1Na1ni0 + Ai1(:,isd,iflow,jflow)*shp1(:,j)*ni0(:,isd)
                  do jsd = 1,nsd
                    Kij1Naj1ni0 = Kij1Naj1ni0 + Kij1(:,isd,jsd,iflow,jflow)*shg1(:,j,jsd)*ni0(:,isd)
                    do kflow = 1,nflow
                      Kij0Nbj0CNa1ni1 = Kij0Nbj0CNa1ni1 + Kij0(:,isd,jsd,iflow,kflow) !AC
     &                                                  * shg0(:,i,jsd)
     &                                                  * cmtrx(:,kflow,jflow)
     &                                                  * shp1(:,j)
     &                                                  * ni1(:,isd)
                    enddo
                  enddo
!E! NOTE: Have to check this (penalty term...)
                  do kflow = 1,nflow
                      CNb0ni0CNa1ni1 = CNb0ni0CNa1ni1 + cmtrx(:,iflow,kflow) * shp0(:,i) * ni0(:,isd) 
     &                                                * cmtrx(:,kflow,jflow) * shp1(:,j) * ni1(:,isd)
                  enddo
                enddo
c
                egmass(:,il,jl) = egmass(:,il,jl)
     &          - (    
     &              pt50 * shp0(:,i) * (Ai1Na1ni0 + Kij1Naj1ni0)
     &          +   pt50 * s         * Kij0Nbj0CNa1ni1
     &          +   e * this_mu(:,iflow,jflow) /length_h * CNb0ni0CNa1ni1
     &            ) * WdetJ0
c
              enddo
            enddo
c
          enddo
        enddo
c
      end subroutine calc_egmass_
c
      subroutine calc_egmass( egmass00, egmass01,
     &                         A0_0, A0_1, Ai0, Ai1,
     &                         Kij0, Kij1,
     &                         AiNa0_, AiNa1_, KijNaj0_, KijNaj1_,
     &                         KijNajC0, KijNajC1,
     &                         shp0, n0, n1, WdetJ0,
     &                         prop,
     &                         nshl0, nshl1)
c        
        real*8, dimension(:,:,:), intent(inout) :: egmass00, egmass01
        real*8, dimension(:,:,:,:), intent(in)  :: AiNa0_, AiNa1_, KijNaj0_, KijNaj1_, KijNajC0, KijNajC1
        real*8, dimension(:,:,:), intent(in):: A0_0,A0_1
        real*8, dimension(:,:,:,:), intent(in) :: Ai0, Ai1
        real*8, dimension(:,:,:,:,:), pointer, intent(in) :: Kij0, Kij1
        real*8, dimension(:,:), intent(in) :: shp0,n0, n1
        real*8, dimension(:), intent(in) :: WdetJ0
        type(prop_t), dimension(:), pointer, intent(in) :: prop
        integer, intent(in) :: nshl0,nshl1
c
        integer :: i,j,p,q,inode0,inode1,isd,jsd
        integer :: i0,j0,il,jl,iflow,jflow
        real*8 :: this_mu(npro,nflow,nflow)
        real*8 :: AiNa0(npro), KijNaj0(npro), KijNbjCn0(npro)
        real*8 :: AiNa1(npro), KijNaj1(npro), KijNbjCn1(npro)
c
c---------> ATTENTION <--------
c  This still needs to be confirmed...
c
      this_mu = zero
      this_mu(:,2,2) = mu(:,2) ! mu
      this_mu(:,3,3) = mu(:,3) ! "
      this_mu(:,4,4) = mu(:,4) ! "
      this_mu(:,5,5) = mu(:,5) ! kappa
c
c---------> END ATTANETION <-----
c
        do i = 1,nshl0
          i0 = nflow*(i-1)
c
          do j = 1,nshl0
            j0 = nflow*(j-1)
c
            do iflow = 1,nflow
              do jflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                AiNa0 = zero
                KijNaj0 = zero
                KijNbjCn0 = zero
c
                do isd = 1,nsd
                  AiNa0 = AiNa0 + pt50 * Ai0(:,isd,iflow,jflow) * n0(:,isd) * shp0(:,j)
                  KijNbjCn0 = KijNbjCn0 + pt50 * s * KijNajC0(:,isd,iflow,jflow) * n0(:,isd)
                  do jsd = 1,nsd
                    KijNaj0 = KijNaj0 + pt50 * Kij0(:,isd,jsd,iflow,jflow) * shg0(:,j,jsd) * n0(:,isd)
                  enddo
                enddo
c
                egmass00(:,il,jl) = egmass00(:,il,jl)
     &          - ( 
     &            + AiNa0 
     &            - KijNaj0
     &            + KijNbjCn0 * shp0(:,j)
CCc     &            + pt50 * A0_0(:,iflow,jflow)*alpha_LF*shp0(:,j)
cc     &            + pt50 * (A0_0(:,iflow,jflow)+A0_1(:,iflow,jflow))*alpha_LF*shp0(:,j)
     &            + e*this_mu(:,iflow,jflow)/length_h * ctc(:,iflow,jflow) * shp0(:,j)
     &            ) * shp0(:,i) * WdetJ0
              enddo
            enddo
c
          enddo
c
          do j = 1,nshl1
            j0 = nflow*(j-1)
c
            do iflow = 1,nflow
              do jflow = 1,nflow
c
                il = i0 + iflow
                jl = j0 + jflow
c
                AiNa1 = zero
                KijNaj1 = zero
                KijNbjCn0 = zero
c
                do isd = 1,nsd
                  AiNa1 = AiNa1 + pt50 * Ai1(:,isd,iflow,jflow) * n0(:,isd) * shp1(:,j)
                  KijNbjCn0 = KijNbjCn0 + pt50 * s * KijNajC0(:,isd,iflow,jflow) * n0(:,isd)
                  do jsd = 1,nsd
                    KijNaj1 = KijNaj1 + pt50 * Kij1(:,isd,jsd,iflow,jflow) * shg1(:,j,jsd) * n0(:,isd)
                  enddo
                enddo
c
                egmass01(:,il,jl) = egmass01(:,il,jl)
     &          - ( 
     &            + AiNa1
     &            - KijNaj1
     &            - KijNbjCn0 * shp1(:,j)
CCc     &            - pt50 * A0_1(:,iflow,jflow)*alpha_LF*shp1(:,j)
cc     &            - pt50 * (A0_0(:,iflow,jflow)+A0_1(:,iflow,jflow))*alpha_LF*shp1(:,j)
     &            - e*this_mu(:,iflow,jflow)/length_h * ctc(:,iflow,jflow) * shp1(:,j)
     &            ) * shp0(:,i) * WdetJ0
              enddo
            enddo
c
          enddo
c
        enddo
c
      end subroutine calc_egmass
c
      end module e3if_lhs_m
