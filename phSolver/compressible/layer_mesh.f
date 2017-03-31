        module layer_mesh_data
          integer layerCommuFlag

        end module
c
        subroutine calc_gc_normal ( x,       shpb,
     &                   ienb,  iBCB,  normal)
c
        include "common.h"
c
        dimension xlb(npro,nenl,nsd),    bnorm(npro,nsd),
     &            rl(npro,nshl,nsd),     WdetJb(npro),
     &            Wfactor(npro)
c
        dimension normal(nshg, nsd)
c
        dimension x(numnp,nsd),
     &            shpb(nshl,ngaussb),      shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),
     &            iBCB(npro,ndiBCB)
c
        dimension lnode(27),               sgn(npro,nshl),
     &            shpfun(npro,nshl),        shdrv(npro,nsd,nshl)
c
        dimension dxdxib(npro,nsd,nsd),      temp(npro),
     &            temp1(npro),               temp2(npro),
     &            temp3(npro),
     &            v1(npro,nsd),              v2(npro,nsd)
c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
c
c.... get the boundary element residuals
c
        rl  = zero
c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d
c     boundary face.
c
c.... set curl to be outward for all types of elememt
c
       if(lcsyst.eq.1) then
         ipt2=2
         ipt3=3
       elseif(lcsyst.eq.2) then
         ipt2=4
         ipt3=2
       elseif(lcsyst.eq.3) then
         ipt2=3
         ipt3=2
       elseif(lcsyst.eq.4) then
         ipt2=2
         ipt3=4
       elseif(lcsyst.eq.5) then
         ipt2=4
         ipt3=2
       elseif(lcsyst.eq.6) then
         ipt2=2
         ipt3=5
       endif
       v1 = xlb(:,ipt2,:) - xlb(:,1,:)
       v2 = xlb(:,ipt3,:) - xlb(:,1,:)
c
c compute cross product
c
       temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
       temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
       temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
c
c mag is area for quads, twice area for tris
c
       temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
       bnorm(:,1) = temp1 * temp
       bnorm(:,2) = temp2 * temp
       bnorm(:,3) = temp3 * temp
c
c
       if (lcsyst .eq. 1) then
         Wfactor(:) = one / (four*temp(:))
       elseif (lcsyst .eq. 2) then
         Wfactor(:) = one / (four*temp(:))
       elseif (lcsyst .eq. 3) then
         Wfactor(:) = one / (two*temp(:))
       elseif (lcsyst .eq. 4) then
         Wfactor(:) = one / (four*temp(:))
       elseif (lcsyst .eq. 5) then
         Wfactor(:) = one / (four*temp(:))
       elseif (lcsyst .eq. 6) then
         Wfactor(:) = one / (two*temp(:))
       endif
c
c.... collect wedge tri and surf ID = BLbaseSrfID
c
       if ((useBLbaseSrfID .eq. 1) .and. (lcsyst .ne. itp_wedge_tri)) then
         do inode = 1,npro
           if (iBCB(inode,2) .ne. BLbaseSrfID) then
             Wfactor(inode) = zero
           endif
         enddo
       endif
c
c.... loop through the integration points
c
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
           shglb=zero  ! protect debugger
           call getshpb(shpb,        shglb,        sgn,
     &                  shpfun,       shdrv)
c
           WdetJb(:) = Qwtb(lcsyst,intp) * Wfactor(:)
c
c  Now lets calculate Integral N_(a:e)^i n_i d Gamma
c
           do n = 1, nshlb
              nodlcl = lnode(n)
              rl(:,nodlcl,1) = rl(:,nodlcl,1) +
     &             shpfun(:,nodlcl) * bnorm(:,1) * WdetJb(:)
              rl(:,nodlcl,2) = rl(:,nodlcl,2) +
     &             shpfun(:,nodlcl) * bnorm(:,2) * WdetJb(:)
              rl(:,nodlcl,3) = rl(:,nodlcl,3) +
     &             shpfun(:,nodlcl) * bnorm(:,3) * WdetJb(:)
           enddo

        enddo  ! quadrature point loop
c
c.... turn outward normal to inward
c
        rl = -rl
c
c.... assemble the normal vector
c
        call local (normal,    rl,     ienb,   3,  'scatter ')
c
c.... end
c
        return
        end

c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine assign_bl_bc ( disp,  iBC, BC, basevID, vID)
c
        include "common.h"
c
c.... please only pass mesh elas BC (:, ndof+2:ndof+5) into this subroutine
c
        dimension disp(numnp,nsd),  iBC(nshg), BC(nshg,4)
c
        integer   basevID, vID
c
            if (btest(iBC(basevID),14)) then
              iBC(vID) = ibset(iBC(vID), 14)
            else
              iBC(vID) = ibclr(iBC(vID), 14)
            endif
c
            if (btest(iBC(basevID),15)) then
              iBC(vID) = ibset(iBC(vID), 15)
            else
              iBC(vID) = ibclr(iBC(vID), 15)
            endif
c
            if (btest(iBC(basevID),16)) then
              iBC(vID) = ibset(iBC(vID), 16)
            else
              iBC(vID) = ibclr(iBC(vID), 16)
            endif
c
            BC(vID, :) = BC(basevID, :)
c
            select case (ibits(iBC(vID),14,3))
            case (1) ! x1 direction
              BC(vID, 1) = disp(vID,1) / Delt(1)
            case (2) ! x2 direction
              BC(vID, 2) = disp(vID,2) / Delt(1)
            case (3) ! x1 & x2 direction
              BC(vID, 1) = disp(vID,1) / Delt(1)
              BC(vID, 3) = disp(vID,2) / Delt(1)
            case (4) ! x3 direction
              BC(vID, 3) = disp(vID,3) / Delt(1)
            case (5) ! x1 & x3 direction
              BC(vID, 1) = disp(vID,1) / Delt(1)
              BC(vID, 3) = disp(vID,3) / Delt(1)
            case (6) ! x2 & x3 direction
              BC(vID, 1) = disp(vID,2) / Delt(1)
              BC(vID, 3) = disp(vID,3) / Delt(1)
            case (7) ! x1 & x2 & x3 direction
              BC(vID, 1) = disp(vID,1) / Delt(1)
              BC(vID, 2) = disp(vID,2) / Delt(1)
              BC(vID, 3) = disp(vID,3) / Delt(1)
            case DEFAULT
              call error('assign_bl_bc','no BC on this point',vID)
            end select
c.... end
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine iBCelas_to_dbl(iBCelas, flag)
c
        real*8 flag
        integer iBCelas
c
        flag = REAL(ibits(iBCelas,14,3))
c
c.... end
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine dbl_to_iBCelas(flag, iBCelas)
c
        real*8 flag
        integer iBCcase, iBCelas
c
        iBCcase = INT(flag + 0.5)
c
            select case (iBCcase)
            case (1) ! x1 direction
              iBCelas = ibset(iBCelas, 14)
              iBCelas = ibclr(iBCelas, 15)
              iBCelas = ibclr(iBCelas, 16)
            case (2) ! x2 direction
              iBCelas = ibclr(iBCelas, 14)
              iBCelas = ibset(iBCelas, 15)
              iBCelas = ibclr(iBCelas, 16)
            case (3) ! x1 & x2 direction
              iBCelas = ibset(iBCelas, 14)
              iBCelas = ibset(iBCelas, 15)
              iBCelas = ibclr(iBCelas, 16)
            case (4) ! x3 direction
              iBCelas = ibclr(iBCelas, 14)
              iBCelas = ibclr(iBCelas, 15)
              iBCelas = ibset(iBCelas, 16)
            case (5) ! x1 & x3 direction
              iBCelas = ibset(iBCelas, 14)
              iBCelas = ibclr(iBCelas, 15)
              iBCelas = ibset(iBCelas, 16)
            case (6) ! x2 & x3 direction
              iBCelas = ibclr(iBCelas, 14)
              iBCelas = ibset(iBCelas, 15)
              iBCelas = ibset(iBCelas, 16)
            case (7) ! x1 & x2 & x3 direction
              iBCelas = ibset(iBCelas, 14)
              iBCelas = ibset(iBCelas, 15)
              iBCelas = ibset(iBCelas, 16)
            case DEFAULT
              call error('dbl_to_iBCelas','wrong BC case',vID)
            end select
c
c
c.... end
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine layerCommuAssembly(global, rtemp, ilwork, n)
c
         include "common.h"
         include "auxmpi.h"
c
         dimension global(nshg,n),
     &             rtemp(maxfront*n,maxtask),
     &             ilwork(nlwork)
c
         real*8  locFlag, revFlag
c
         integer itkbeg, jdl, iacc, numseg, isgbeg, icsg, icid,
     &           lenseg, itemp, lfront, j, idof
c
         numtask = ilwork(1)
         itkbeg=1 ! slot in ilwork
         jdl=0    ! slot in rtemp
c
         do j=1,numtask        ! loop over all tasks
c
c.... total number of nodes involved in this task (lfront)
c
            lfront = 0
            do is = 1,numseg
              lenseg = ilwork (itkbeg + 4 + 2*is)
              lfront = lfront + lenseg
            enddo
c
            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)
            isgbeg = ilwork (itkbeg + 5)
            if(iacc.eq.1) then
               jdl=jdl+1  ! keep track of order of rtemp's
c
c.... only add the data from growth curve once
c
               itemp = 1
               do idof = 1,3
                 do is = 1,numseg
                   isgbeg = ilwork (itkbeg + 3 + 2*is)
                   lenseg = ilwork (itkbeg + 4 + 2*is)
c                 isgend = isgbeg + lenseg - 1
                   do icsg = 1,lenseg
                     icid  = isgbeg + icsg - 1
                     locFlag = global(icid,4)
                     revFlag = rtemp(3*lfront+mod(itemp,lfront),jdl)
                     if(revFlag .gt. 0.5 .and. locFlag .lt. 0.5) then
                       global(icid,idof) = global(icid,idof)
     &                                   + rtemp (itemp,jdl)
                       global(icid,4) = revFlag
                     endif
                     itemp = itemp + 1
                   enddo
                 enddo
               enddo
            endif ! end of receive (iacc=1)
            itkbeg = itkbeg + 4 + 2*numseg
         enddo
c
c.... end
c
        return
        end
c

