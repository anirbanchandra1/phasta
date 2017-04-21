        subroutine ElmGMRElas (x,       disp,    shp,     shgl,
     &                         iBC,     BC,      shpb,    shglb,
     &                         shpif,   elasres, elasBDiag,
     &                         iper,    ilwork,  elaslhsK,
     &                         col,     row,     meshq)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual
c vector, and the preconditioning matrix, for mesh-elastic solve
c
c----------------------------------------------------------------------
c
         use pointer_data
         use timedataC
         use readarrays ! read BLflt, BLgr, BLtnv, BLlist
c
        include "common.h"
        include "mpif.h"
c
        integer col(nshg+1), row(nnz*nshg)
        real*8  elaslhsK(nelas*nelas,nnz_tot),
     &          meshq(numel),
     &          meshV(numel)
c
        real*8  gcnormal(nshg, nsd)
c
        dimension x(numnp,nsd),        disp(numnp,nsd),
     &            xtmp(numnp,nsd),     iBC(nshg),
     &            BC(nshg,ndofBC),
     &            elasres(nshg,nelas),
     &            elasBDiag(nshg,nelas,nelas),
     &            iper(nshg)
c
        real*8  flag(numnp),   dispCommu(numnp,4)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)
c
        real*8, dimension(maxtop,    maxsh,maxqpt) :: shpif
c
        dimension ilwork(nlwork)
c
        integer errorcount(2)
c
        integer listcounter, ngc, itnv, basevID, nv, vID, vID2, ioffset
        real*8  iflt, igr, igrexp, tmp
        real*8  inormal(nsd)
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: Estiff(:,:,:)
c
c.... ------------------->   check layered mesh   <--------------------
c
        real*8, allocatable :: tmpNm(:,:)
        real*8, allocatable :: thickness(:), gcvpos(:)
        real*8  tmpCP1,tmpCP2,tmpCP3,norm1,norm2
c
c        do i = 1,numnp
c          if( x(i,1) .ge. 0.0043890 .and. x(i,1) .le. 0.0043892) then
c          if( x(i,1) .ge. -0.021889 .and. x(i,1) .le. -0.021887) then
c          if( x(i,1) .ge. -2.40772e-2 .and. x(i,1) .le. -2.40770e-2) then
c            write(*,*) "find base: rank:",myrank,"ID",i
c            write(*,*) "find top most: rank:",myrank,"ID",i
c          endif
c        enddo
c
c... check top most
         if (myrank .eq. 0) then
           write(*,*) "top most rank:",myrank,"coord:",x(1466,:)
         endif
         if (myrank .eq. 2) then
           write(*,*) "top most rank:",myrank,"coord:",x(2729,:)
         endif
         if (myrank .eq. 3) then
           write(*,*) "top most rank:",myrank,"coord:",x(2757,:)
         endif
c... not end check top most
c
c.... loop over growth curves
c
c        listcounter = 0
c        ioffset = 1 ! the ID starts from 1 in phasta
c        do ngc = 1, numgc
c         itnv = BLtnv(ngc) ! number of vertices on this growth curve
c          allocate(tmpNm(itnv-1, nsd))
c          allocate(thickness(itnv-1))
c          allocate(gcvpos(itnv-1))
c          basevID = BLlist(listcounter + 1) + ioffset
c          do nv = 2, itnv
c            vID = BLlist(listcounter + nv) + ioffset
c
c... continue check top most
c            if(myrank .eq. 2 .and. basevID .eq. 393 .and. nv .eq. 4) then
c              write(*,*) "top most: ",x(vID,:)
c            endif
c... end check top most
c
c... check inter part gc positions
c            if((myrank .eq. 0 .and. basevID .eq. 555) .or.
c     &         (myrank .eq. 2 .and. basevID .eq. 14 ) .or.
c     &         (myrank .eq. 3 .and. basevID .eq. 25 )) then
c              gcvpos(nv-1) = x(vID, 3) ! 1,2,3
c            endif
c... not end check inter part gc positions
c
c            vID2= BLlist(listcounter + nv - 1) + ioffset ! the previous one
c            tmpNm(nv-1,:) = x(vID,:) - x(vID2,:)
c
c... check thickness
c            thickness(nv-1) = sqrt(tmpNm(nv-1,1)*tmpNm(nv-1,1)
c     &                           + tmpNm(nv-1,2)*tmpNm(nv-1,2)
c     &                           + tmpNm(nv-1,3)*tmpNm(nv-1,3))
c
c... not end check thickness
c
c
c... check normal direction alignment
c            if(nv .gt. 2) then
c              tmpCP1 = tmpNm(nv-2,2)*tmpNm(nv-1,3)-tmpNm(nv-1,2)*tmpNm(nv-2,3)
c              tmpCP2 = tmpNm(nv-1,1)*tmpNm(nv-2,3)-tmpNm(nv-2,1)*tmpNm(nv-1,3)
c              tmpCP3 = tmpNm(nv-2,1)*tmpNm(nv-1,2)-tmpNm(nv-1,1)*tmpNm(nv-2,2)
c              norm1  = sqrt(tmpNm(nv-2,1)*tmpNm(nv-2,1)
c     &                    + tmpNm(nv-2,2)*tmpNm(nv-2,2)
c     &                    + tmpNm(nv-2,3)*tmpNm(nv-2,3))
c              norm2  = sqrt(tmpNm(nv-1,1)*tmpNm(nv-1,1)
c     &                    + tmpNm(nv-1,2)*tmpNm(nv-1,2)
c     &                    + tmpNm(nv-1,3)*tmpNm(nv-1,3))
c              write(*,*) "rank: ",myrank," gc: ",ngc," layer: ",nv-1,
c     &          " cp/m: ",tmpCP1/norm1/norm2,tmpCP2/norm1/norm2,tmpCP3/norm1/norm2
c            endif
c... end check normal direction alignment
c
c          enddo ! over this growth curve
c
c... continue check inter part gc positions
c          if((myrank .eq. 0 .and. basevID .eq. 555) .or.
c     &         (myrank .eq. 2 .and. basevID .eq. 14 ) .or.
c     &         (myrank .eq. 3 .and. basevID .eq. 25 )) then
c            write(*,*) "rank:",myrank,"gc:",ngc,"v x pos:",gcvpos(:)
c          endif
c... end check inter part gc positions
c
c... continue check thickness
c          write(*,*) "rank: ",myrank," gc: ",ngc," thickness: ",thickness(:)
c... end check thickness
c
c          deallocate(gcvpos)
c          deallocate(thickness)
c          deallocate(tmpNm)
c          listcounter = listcounter + itnv ! update counter
c        enddo
c
c.... ------------------->   layer base elements   <-------------------
c
c.... calculate the normal of each growth curve based on new boundary positions
c
        xtmp = x + disp
        gcnormal = zero
c
c.... loop over the boundary elements
c
        boundary_blocks: do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel
c
          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... collect wedge_tri or surfID option is on
c
          if((lcsyst.ne.itp_wedge_tri) .and. (useBLbaseSrfID.eq.0)) cycle
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,             shpb(lcsyst,1:nshl,:),
     &                         mienb(iblk)%p,    miBCB(iblk)%p,
     &                         gcnormal)
c
        enddo boundary_blocks ! end loop the boundary elements
c
c.... loop over the interface elements
c
        interface_blocks: do iblk = 1, nelblif
c
          iel     = lcblkif(1, iblk)
          npro    = lcblkif(1,iblk+1) - iel
          lcsyst0 = lcblkif(3, iblk)    ! element0 type
          lcsyst1 = lcblkif(4, iblk)    ! element1 type
          ipord   = lcblkif(5, iblk)    ! polynomial order
          nenl0   = lcblkif(6, iblk)    ! number of vertices per element0
          nenl1   = lcblkif(7, iblk)    ! number of vertices per element1
          mater0  = lcblkif(9, iblk)
          mater1  = lcblkif(10,iblk)
          nshl0   = lcblkif(iblkif_nshl0,iblk)
          nshl1   = lcblkif(iblkif_nshl1,iblk)
          itpid   = lcblkif(iblkif_topology,iblk)
          ngaussif = nintif(itpid)
c
c.... the 0 side
c
          lcsyst = lcsyst0
          nenl = nenl0
          nshl = nshl0

c          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... collect wedge_tri or surfID option is on
c
          if((lcsyst.ne.itp_wedge_tri) .and. (useBLbaseSrfID.eq.0)) cycle
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,            shpif(lcsyst,1:nshl,:),
     &                         mienif0(iblk)%p, miBCB(iblk)%p,
     &                         gcnormal)
c
c.... end of the 0 side
c
c
c.... the 1 side
c
          lcsyst = lcsyst1
          nenl = nenl1
          nshl = nshl1

c          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... collect wedge_tri or surfID option is on
c
          if((lcsyst.ne.itp_wedge_tri) .and. (useBLbaseSrfID.eq.0)) cycle
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,            shpif(lcsyst,1:nshl,:),
     &                         mienif1(iblk)%p, miBCB(iblk)%p,
     &                         gcnormal)
c
c.... end of the 1 side
c
        enddo interface_blocks ! end loop the interface elements
c
c.... communication
c
        if (numpe > 1) then
          call commu (gcnormal  , ilwork, nsd  , 'in ')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          call commu (gcnormal  , ilwork, nsd  , 'out ')
        endif
c
c.... end calculation of growth curve normal
c
c.... ---------------->   Re-position layered mesh   <-----------------
c
c.... loop over growth curves
c
        flag = zero
        listcounter = 0
        ioffset = 1 ! the ID starts from 1 in phasta
        do ngc = 1, numgc
          itnv = BLtnv(ngc) ! number of vertices on this growth curve
c
c.... precaution
c
          if (itnv .lt. 2) then
            listcounter = listcounter + itnv
            cycle ! not loop over vertices
          endif
c
c.... prepare other paramteres
c
          iflt = BLflt(ngc) ! first layer thickness of this growth curve
          igr  = BLgr(ngc)  ! growth ratio of this growth curve
          basevID = BLlist(listcounter + 1) + ioffset
          tmp  = sqrt( gcnormal(basevID,1) * gcnormal(basevID,1)
     &               + gcnormal(basevID,2) * gcnormal(basevID,2)
     &               + gcnormal(basevID,3) * gcnormal(basevID,3) )
          inormal(:) = gcnormal(basevID,:) / tmp
c
c.... loop over vertices on this growth curve
c
          igrexp = 1.0
          do nv = 2, itnv
            vID = BLlist(listcounter + nv) + ioffset
            vID2= BLlist(listcounter + nv - 1) + ioffset ! the previous one
            xtmp(vID,:) = xtmp(vID2,:) + iflt * inormal(:) * igrexp  ! igr**(nv-2)
            disp(vID,:) = xtmp(vID,:) - x(vID,:)
            igrexp = igrexp * igr
            flag(vID) = 1.0 ! setup flag
          enddo ! over this growth curve
c
          listcounter = listcounter + itnv ! update counter
c
c.... end loop growth curves
c
        enddo
c
c.... end re-position layered mesh
c
c.... -------------------->     communication     <--------------------
c
c.... combine disp + flag into dispCommu
c
        do i = 1, numnp
          dispCommu(i,1:3) = disp(i,:)
          dispCommu(i,4)   = flag(i)
        enddo
c
c.... turn on commu flag
c
        layerCommuFlag = 1
c
c.... commu
c
        if (numpe > 1) then
          call commu (dispCommu  , ilwork, 4  , 'in ')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          call commu (dispCommu  , ilwork, 4  , 'out ')
        endif
c
c.... turn off commu flag
c
        layerCommuFlag = 0
c
c.... separate dispCommu into disp; update BC
c
        do i = 1, numnp
          disp(i,:) = dispCommu(i,1:3)
          call setBLbc(disp(i,:), iBC(i), BC(i, ndof+2:ndof+5))
        enddo
c
c.... -------------------->   interior elements   <--------------------
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        elasres = zero
        errorcount = zero
        if (lhs. eq. 1)    elaslhsK  = zero
        if (iprec .ne. 0)  elasBDiag = zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk            ! used in timeseries
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
          ndofelas = nshl * nelas
c
          allocate (Estiff(npro,ndofelas,ndofelas))
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
c
          Estiff = zero
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c
c.... Shape measure. Calculate the shape quality
c
          call shpMeasure(x, mien(iblk)%p, tmpshp, tmpshgl,
     &                    meshq(iel:iel+npro-1),
     &                    meshV(iel:iel+npro-1), errorcount )
c
          call AsIGMRElas (x,             disp,
     &                     tmpshp,        tmpshgl,
     &                     mien(iblk)%p,  elasres,
     &                     elasBDiag,     Estiff,
     &                     meshq(iel:iel+npro-1),
     &                     meshV(iel:iel+npro-1)   )
c
c.... satisfy the BCs on the implicit LHS
c
          call bc3LHSElas (iBC, BC(:, ndof+2:ndof+5),
     &                     mien(iblk)%p, Estiff)
c
c.... Fill-up the global sparse LHS mass matrix
c
          call fillsparseElas( mien(iblk)%p, Estiff,
     &                         elaslhsK, row, col)
c
          deallocate ( Estiff )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
c
      if (errorcount(1) .lt. 0 .or. errorcount(2) .lt. 0) then
        write(*,*) errorcount(1), " elements Meshq larger than one; ",
     &             errorcount(2), " elements Meshq smaller than zero."
      endif
c
      if(iabc==1) then               ! are there any axisym BCs
          call rotabc(elasres(1,2), iBC,  'in ')
      endif
c
c.... -------------------->   communications <-------------------------
c
      if (numpe > 1) then
          call commu (elasres  , ilwork, nelas  , 'in ')
c
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c
          if (iprec .ne. 0) then
             call commu (elasBDiag,    ilwork,
     &                   nelas*nelas,  'in ')
          endif
      endif
c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResElas (iBC,     BC(:, ndof+2:ndof+5),
     &                 elasres, iper,    ilwork)
c
c.... satisfy the BCs on the block-diagonal preconditioner
c
      if (iprec .ne. 0) then
         call bc3BDgElas (iBC,       BC(:, ndof+2:ndof+5),
     &                    elasBDiag, iper,    ilwork)
      endif
c
c.... return
c
      return
      end

