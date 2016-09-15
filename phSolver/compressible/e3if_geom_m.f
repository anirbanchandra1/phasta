      module e3if_geom_m
c
        use e3if_defs_m
        use hierarchic_m
c
        implicit none
c
      contains
c
        subroutine malloc_e3if_geom
c
          allocate(xl0(npro,nenl0,nsd))
          allocate(xl1(npro,nenl1,nsd))
          allocate(area(npro))
          allocate(WdetJif0(npro),WdetJif1(npro))
          allocate(if_normal_l0(npro,nshl0,nsd+1))
          allocate(if_normal_l1(npro,nshl1,nsd+1))
          allocate(nv0(npro,nsd))
          allocate(nv1(npro,nsd))
          allocate(shp0(npro,nshl0))
          allocate(shp1(npro,nshl1))
          allocate(shgl0(npro,nsd,nshl0))
          allocate(shgl1(npro,nsd,nshl1))
          allocate(sgn0(npro,nshl0))
          allocate(sgn1(npro,nshl1))
c
        end subroutine malloc_e3if_geom
c
        subroutine mfree_e3if_geom
c
          deallocate(xl0,xl1)
          deallocate(area)
          deallocate(WdetJif0,WdetJif1)
          deallocate(if_normal_l0,if_normal_l1)
          deallocate(nv0,nv1)
          deallocate(shp0,shp1)
          deallocate(shgl0,shgl1)
          deallocate(sgn0,sgn1)
c
        end subroutine mfree_e3if_geom
c
        subroutine calc_if_normals(xl0,xl1,shpif0,shpif1,shgif0,shgif1,qwtif0,qwtif1)
c
          real*8,  dimension(:,:,:), pointer, intent(in) :: xl0, xl1 
          real*8, dimension(nshl0,nqpt),intent(in)   :: shpif0
          real*8, dimension(nshl1,nqpt),intent(in)   :: shpif1
          real*8, dimension(nsd,nshl0,nqpt), intent(in) :: shgif0
          real*8, dimension(nsd,nshl1,nqpt), intent(in) :: shgif1
          real*8, dimension(nqpt), intent(in) :: qwtif0, qwtif1
c
          integer :: intp
c
          do intp = 1, nqpt
c
c... calculate normal vectors at the integration point...
c
            call calc_normal_vectors(nv0,area,WdetJif0,xl0,qwtif0,lcsyst0,intp,npro)
            call calc_normal_vectors(nv1,area,WdetJif1,xl1,qwtif1,lcsyst1,intp,npro)
c
            call  getshp(shp0, shgl0, shpif0, shgif0, sgn0, npro, nsd, nshl0, nqpt, nenl0, intp, ipord)
            call  getshp(shp1, shgl1, shpif1, shgif1, sgn1, npro, nsd, nshl1, nqpt, nenl1, intp, ipord)
c
c... distribute the weighted vectors from integration points to nodes
c
            call collect_weighted_normals(if_normal_l0,nv0,shp0,WdetJif0,nshl0)
            call collect_weighted_normals(if_normal_l1,nv1,shp1,WdetJif1,nshl1)
c
          enddo
c
        end subroutine calc_if_normals
c
        subroutine collect_weighted_normals(if_normal_l,nv,shp,WdetJif,nshl)
c
          real*8, dimension(:,:,:), intent(inout) :: if_normal_l
          real*8, dimension(:,:), intent(in) :: nv
          real*8, dimension(:,:),   intent(in)    :: shp
          real*8, dimension(:),     intent(in)    :: WdetJif
          integer, intent(in) :: nshl
c
          integer :: n,isd
c
          do n = 1,nshl
            do isd = 1,nsd
              if_normal_l(:,n,isd) = if_normal_l(:,n,isd) 
     &                               + shp(:,n)*nv(:,isd)*area(:)
c     &                               + shp(:,n)*nv(:,isd)*WdetJif(:)
            enddo
          enddo
c
        end subroutine collect_weighted_normals
c
        subroutine calc_normal_vectors(nv,area,WdetJ,xl,qwt,lcsyst,intp,npro)
c
          real*8, dimension(:,:), pointer, intent(out) :: nv
          real*8, dimension(:),   pointer, intent(out) :: area, WdetJ
          real*8, dimension(:,:,:), pointer, intent(in) :: xl
          real*8, dimension(nqpt),           intent(in) :: qwt
          integer, intent(in) :: lcsyst,intp,npro
c
          real*8 :: temp_len(npro)
          real*8, dimension(npro, nsd) :: v1, v2, temp_normal
          integer :: isd, iel
          character(len=8) :: err_msg
c
c      write(*,*) 'In calc_normal_vectors...'
c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
          do isd = 1,nsd
            v1(:,isd) = xl(:,2,isd) - xl(:,1,isd)
            v2(:,isd) = xl(:,3,isd) - xl(:,1,isd)
          enddo
c
          select case (lcsyst)
          case (1)
c
            temp_normal(:,1) = + v1(:,2)*v2(:,3) - v1(:,3)*v2(:,2)
            temp_normal(:,2) = - v1(:,1)*v2(:,3) + v1(:,3)*v2(:,1)
            temp_normal(:,3) = + v1(:,1)*v2(:,2) - v1(:,2)*v2(:,1)
c
            area = pt50 * sqrt(temp_normal(:,1)*temp_normal(:,1)
     &                       + temp_normal(:,2)*temp_normal(:,2)
     &                       + temp_normal(:,3)*temp_normal(:,3))
c
          case default
c            write(err_msg,'(a,i1)') 'lcsyst ',lcsyst
            call error ('calc_normal_vectors ', err_msg, 0)
          end select
c
          do iel = 1,npro
            temp_len(iel) = sqrt(dot_product(temp_normal(iel,1:3),temp_normal(iel,1:3)))
            nv(iel,1:nsd) = temp_normal(iel,1:nsd) / temp_len(iel)
          enddo
c
c... also calculate WdetJ here...
c
          select case (lcsyst)
          case (1)
            WdetJ = qwt(intp) * temp_len * pt25
          case default
            write(err_msg,'(a,i1)') 'lcsyst ',lcsyst
            call error ('calc WdetJif ', err_msg, 0)
          end select
c
        end subroutine calc_normal_vectors
c
        subroutine calc_mean_curvature(if_kappa,ienif0,ienif1)
c
          real*8, dimension(:,:), pointer, intent(inout) :: if_kappa
          integer, dimension(:,:), pointer, intent(in) :: ienif0,ienif1
c
          integer,parameter :: nnode = 3 ! ONLY WORKS ON THE TRIAGLES
          integer :: inode(3,3)
          data inode  /1,2,3,2,3,1,3,1,2/
          integer :: iel,i,j,k
          real*8 :: v1(3),v2(3),v(3),area
          real*8, dimension(nnode) :: cot0(3),cot1(3)
          logical :: obtuse
c
          do iel = 1,npro
c
            do k = 1,nnode
              i = inode(k,2)
              j = inode(k,3)
              v1 = xl0(iel,i,:) - xl0(iel,k,:)
              v2 = xl0(iel,j,:) - xl0(iel,k,:)
              cot0 = cotan(v1,v2)
              v1 = xl1(iel,i,:) - xl1(iel,k,:)
              v2 = xl1(iel,j,:) - xl1(iel,k,:)
              cot1 = cotan(v1,v2)
            enddo
c
            call cross(v,v1,v2)
            area = pt50*norm2(v)
c
            obtuse = any(cot0(:) < 0)
c
            do k = 1,nnode
c
              i = inode(k,2)
              j = inode(k,3)
c
              if_kappa(ienif0(iel,i),1:nsd) = if_kappa(ienif0(iel,i),1:nsd) + pt5*cot0(k)*(xl0(iel,i,1:nsd)-xl0(iel,j,1:nsd))
              if_kappa(ienif0(iel,j),1:nsd) = if_kappa(ienif0(iel,j),1:nsd) + pt5*cot0(k)*(xl0(iel,j,1:nsd)-xl0(iel,i,1:nsd))
c
              if_kappa(ienif1(iel,i),1:nsd) = if_kappa(ienif1(iel,i),1:nsd) + pt5*cot1(k)*(xl1(iel,i,1:nsd)-xl1(iel,j,1:nsd))
              if_kappa(ienif1(iel,j),1:nsd) = if_kappa(ienif1(iel,j),1:nsd) + pt5*cot1(k)*(xl1(iel,j,1:nsd)-xl1(iel,i,1:nsd))
c
c ...NOW collect the voronoi area...
c
              if (.not. obtuse) then
c
                v = xl0(iel,i,:) - xl0(iel,j,:)
                if_kappa(ienif0(iel,k),nsd+1) = if_kappa(ienif0(iel,k),nsd+1) + pt125*cot0(k)*norm2(v)
c
                v = xl1(iel,i,:) - xl1(iel,j,:)
                if_kappa(ienif1(iel,k),nsd+1) = if_kappa(ienif1(iel,k),nsd+1) + pt125*cot1(k)*norm2(v)
c
              else
c
                if (cot0(k) < zero) then
                  if_kappa(ienif0(iel,k),nsd+1) = if_kappa(ienif0(iel,k),nsd+1) + pt50*area
                else
                  if_kappa(ienif0(iel,k),nsd+1) = if_kappa(ienif0(iel,k),nsd+1) + pt25*area
                endif
                if (cot1(k) < zero) then
                  if_kappa(ienif1(iel,k),nsd+1) = if_kappa(ienif1(iel,k),nsd+1) + pt50*area
                else
                  if_kappa(ienif1(iel,k),nsd+1) = if_kappa(ienif1(iel,k),nsd+1) + pt25*area
                endif
              endif
c
            enddo
c
          enddo
c
        end subroutine calc_mean_curvature
c
        subroutine cross(v,v1,v2)
          real*8, dimension(nsd), intent(out) :: v
          real*8, dimension(nsd), intent(in) :: v1,v2
          v(1) = v1(2)*v2(3) - v1(3)*v2(2)
          v(2) = v1(3)*v2(1) - v1(1)*v2(3)
          v(3) = v1(1)*v2(2) - v1(2)*v2(1)
        end subroutine cross
c
        real*8 function cotan(v1,v2)
          real*8, dimension(nsd) :: v1,v2,v
          call cross(v,v1,v2)
          cotan = dot_product(v1,v2)/norm2(v)
        end function cotan
c
      end module e3if_geom_m
