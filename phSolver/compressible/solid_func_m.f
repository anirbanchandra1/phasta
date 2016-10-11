      module solid_func_m
c
        use e3_def_m
        use solid_m
        use eqn_state_m
c
        implicit none
c
        integer :: almBi_s, alfBi_s, gamBi_s
        
c
        real*8  :: Delt_s
c
        real*8, dimension(:,:), pointer :: dudx, dudy, dudz
        real*8, dimension(:,:,:), pointer :: AS
c
      contains
c
      subroutine calc_solid
c
        allocate(AS(npro,6,6))
c
        call calc_as_matrix
        d = almBi * bq + alfBi * Delt(1) * (almBi - gamBi) * bq_dot
        call setB_af(d, AS, bq_af)
        call get_det(bq_af,det_baf)
        Ja_def= (det_baf)**0.5
        call get_det
c
        deallocate(AS)

c
      end subroutine calc_solid
c
      subroutine calc_as_matrix(dudx, dudy,dudz)
c initialize the AS matrix
         AS = zero !add
c
         AS(:,1,1) = 2 * dudx(:,1)
         AS(:,1,5) = 2 * dudz(:,1)
         AS(:,1,6) = 2 * dudy(:,1)
c
         AS(:,2,2) = 2 * dudy(:,2)
         AS(:,2,4) = 2 * dudz(:,2)
         AS(:,2,6) = 2 * dudx(:,2)
c
         AS(:,3,3) = 2 * dudz(:,3)
         AS(:,3,4) = 2 * dudy(:,3)
         AS(:,3,5) = 2 * dudx(:,3)
c
         AS(:,4,2) = dudy(:,3)
         AS(:,4,3) = dudz(:,2)
         AS(:,4,4) = dudy(:,2) + dudz(:,3)
         AS(:,4,5) = dudx(:,2)
         AS(:,4,6) = dudx(:,3)
c
         AS(:,5,1) = dudx(:,3)
         AS(:,5,3) = dudz(:,1)
         AS(:,5,4) = dudy(:,1)
         AS(:,5,5) = dudx(:,1) + dudz(:,3)
         AS(:,5,6) = dudy(:,3)
c         
         AS(:,6,1) = dudx(:,2)
         AS(:,6,2) = dudy(:,1)
         AS(:,6,4) = dudz(:,1)
         AS(:,6,5) = dudz(:,2)
         AS(:,6,6) = dudx(:,1) + dudy(:,2)
c
c
      end subroutine calc_as_matrix
c
      subroutine setB_af (d, AS, bq_af)
c... calculate the left Cauchy-green tensor at time step n+af
       implicit none 
c
       real*8, dimension(npro,6),intent(in) :: d
       real*8, dimension(npro,6,6),intent(in) :: AS
       integer,parameter :: nsize = 6 
c
       real*8, dimension(npro,6) :: bq_af
       real*8, dimension(6,6) :: ident
       real*8, dimension(6,6) :: temp_matrix
c..................
       ident = zero
        do i = 1, nsize
        ident(i,i) = one
        enddo
c
        do i = 1, npro
        temp_matrix(:,:) = almBi * ident(:,:) + gamBi * Delt(1)
     &                     *alfBi * AS(i,:,:) !check here
        bq_af(i,:) = matmul(temp_matrix(:,:) , d(i,:))
        enddo
c
c
       end subroutine setB_af 
c
c
c
c
       subroutine get_det(matr,det_matr)
c... calculate the determinant of 3 by 3 matrix
      implicit none 
c
      real*8, dimension(npro,6),intent(in) :: matr
c
      real*8, dimension(npro) :: det_matr
c..................
      det_matr(:) = matr(:,1) * matr(:,2) * matr(:,3)
     &+             matr(:,6) * matr(:,4) * matr(:,5)
     &+             matr(:,5) * matr(:,6) * matr(:,4) 
     &-             matr(:,1) * matr(:,4) * matr(:,4)
     &-             matr(:,5) * matr(:,2) * matr(:,5)
     &-             matr(:,6) * matr(:,6) * matr(:,3)
c
       end subroutine get_det
c
c
      end module solid_func_m
