      module global_const_m
c
        use iso_c_binding
c
        implicit none
c
c.... parameters  IF YOU CHANGE THES YOU HAVE TO CHANGE THEM IN
c                  common_c.h ALSO
c
        integer, parameter     :: MAXBLK = 50000
        integer, parameter     :: MAXSH = 32, NSD = 3 , NSDSQ = 9
c
c  The five types of region topology are  1= Tet, 2=Hex, 3= Wedge (tri-start),
c                                         4= Wedge (quad-first) 5=pyramid
c
c  The two types of face topology are  1= tri, 2=quad
c
        integer, parameter     :: MAXTOP = 6, MAXSURF=1000
c
c
c...  The twelve different topological interface region are:
c
        integer, parameter     :: MAXTOPIF = 12
c
c  sharing a tri face:
c
c  1= tet-tet 
c  2= tet-pyramid
c  3= tet-wedge
c  4= pyramid-pyramid
c  5= pyramid-wedge
c  6= wedge-wedge
c
c  sharing a quad face:
c
c  7= pyramid-pyramid
c  8= pyramid-wedge
c  9= pyramid-hex
c  10= wedge-wedge
c  11= wedge-hex
c  12= hex-hex
c
c
      integer, parameter :: MAXQPT = 125
c
      end module global_const_m
c
c----------------------------------------------------------------------
c
c.... common /conpar/   : input constants
c
c numnp         : number of nodal points
c numel         : number of elements
c numelb        : number of boundary elements
c numpbc        : number of nodes having a boundary condition
c nen           : maximum number of element nodes
c nfaces        : maximum number of element faces
c nsd           : number of space dimensions
c numflx        : number of flux boundary nodes
c ndof          : number of degrees of freedom per node
c iALE          : ALE formulation flag
c iSOLID        : Solid formulation flag
c icoord        : coordinate system flag
c navier        : Navier-Stokes calculation flag
c irs           : restart option 
c iexec         : execute flag
c necho         : input echo parameter
c ichem         : equilibrium chemistry flag (for outchem.step dump)
c iRK           : Runge-Kutta flag
c nshg          : global number of shape functions (degrees of freedom,
c                 or equations). Computed from the specified p-order,
c                 the number of edges, and the number of faces (in the
c                 entire mesh)
c
c----------------------------------------------------------------------
c
      module conpar_m
c
        use iso_c_binding
c
        implicit none
c
        integer, target :: numnp, numel,  numelb, numelif,
     &                  numpbc,   nen,    nfaces,
     &                  numflx,   ndof,   iALE, iSOLID,
     &                  icoord,   navier,
     &                  irs,      iexec,  necho,  ichem,  iRK,    nedof,
     &                  ndofelas, nshg,   nnz,    istop,  nflow,  nelas, 
     &                  nnz_tot,  idtn,
     &                  ncorpsize, iownnodes, usingpetsc
c
        common /conpar/ numnp, numel,  numelb, numelif,
     &                  numpbc,   nen,    nfaces,
     &                  numflx,   ndof,   iALE, iSOLID,
     &                  icoord,   navier,
     &                  irs,      iexec,  necho,  ichem,  iRK,    nedof,
     &                  ndofelas, nshg,   nnz,    istop,  nflow,  nelas, 
     &                  nnz_tot,  idtn,
     &                  ncorpsize, iownnodes, usingpetsc
c
      end module conpar_m
c
c----------------------------------------------------------------------
c
c.... common /timdat/   : time data
c
c time          : current run time
c CFLfld        : CFL number for fluid flow
c CFLsld        : CFL number for structural heating
c Dtgl          : inverse of global time step
c Dtmax         : maximum delta-time
c alpha         : trapezoidal rule parameter
c etol          : epsilon tolerance for GMRES
c lstep         : current time step
c ifunc         : func. eval. counter (=niter*(lstep-lstep0) + iter)
c itseq         : sequence number
c istep         : step number (reseted at the beginning of the run)
c iter          : iteration number
c nitr          : number of multi-corrector iterations for this sequence
c
      module time_m
c
        use iso_c_binding
c
        implicit none
c
        integer :: lstep, ifunc, itseq, istep, iter, nitr, iCFLworst, lskeep
        real*8 :: time, CFLfld, CFLsld, Dtgl, Dtmax, alpha, etol,
     &            almi, alfi, gami, almBi, alfBi, gamBi, flmpl, flmpr, dtol(2)
c
        common /timdat/ time,   CFLfld, CFLsld, Dtgl,   Dtmax,  alpha,
     &                  etol,   lstep,  ifunc,  itseq,  istep,  iter,
     &                  nitr,   almi,   alfi,   gami,   
     &                  almBi,  alfBi,  gamBi,
     &                  flmpl,  flmpr,
     &                  dtol, iCFLworst, lskeep
c
      end module time_m
c
c----------------------------------------------------------------------
c
c.... common /elmpar/   : element parameters
c
c lelCat        : element category (P1, Q1, P2, Q2, etc.)
c lcsyst        : element coordinate system
c iorder        : element order (=k for Pk and Qk)
c nenb          : number of element nodes per boundary sides
c maxsh         : total number integration points
c maxshb        : total number integration points of boundary elements
c nelblk        : number of element blocks
c nelblb        : number of boundary element blocks
c nelblif       : number of interface element blocks
c ndofl         : number of degrees of freedom (for current block)
c nsymdl        : number of d.o.f for symm. storage (for current block)
c nenl          : number of element nodes (for current block)
c nfacel        : number of element faces (for current block)
c nenbl         : number of boundary element nodes
c intind        : integration data index
c nintg         : number of integration points
c mattyp        : material type ( = 0 for fluid; = 1 for solid )
c iftpid(MAXBLK): holds the interface topological combination
c
      module elmpar_m
c
        use iso_c_binding
        use global_const_m
c
        implicit none
c
	      integer, target ::  nelblk, nelblb, nelblif
        integer :: lelCat, lcsyst, iorder, nenb,   
     &                  ndofl,  nsymdl, nenl,   nfacel,
     &                  nenl0,  nenl1,  lcsyst0, lcsyst1,
     &                  nenbl,  intind, mattyp,
     &                  mattyp0, mattyp1, 
     &                  iftpid(MAXBLK)
c
        common /elmpar/ lelCat, lcsyst, iorder, nenb,   
     &                  nelblk, nelblb, nelblif,
     &                  ndofl,  nsymdl, nenl,   nfacel,
     &                  nenbl,  intind, mattyp,
     &                  iftpid
c
      end module elmpar_m
c
c----------------------------------------------------------------------
c
c.... common /blkdat/   : blocking data
c
c lcblk  (10,MAXBLK+1) : blocking data for the interior elements
c lcblkb (10,MAXBLK+1) : blocking data for the boundary elements
c lcblkif (14,MAXBLK+1) : blocking data for the interface elements
c
      module blkdat_m
c
        use global_const_m
        implicit none
c
        integer :: lcblk  (10,MAXBLK+1),
     &             lcblkb (10,MAXBLK+1),
     &             lcblkif(14,MAXBLK+1)
c
      end module blkdat_m
c
c----------------------------------------------------------------------
c
      module intpt_m
c
        use global_const_m
        implicit none
c
c.... hierarchic basis functions
c
        real*8 :: Qpt (MAXTOP ,4,MAXQPT), Qwt (MAXTOP ,MAXQPT), 
     &            Qptb(MAXTOP,4,MAXQPT),  Qwtb(MAXTOP,MAXQPT), 
     &            Qptif0(MAXTOPIF,4,MAXQPT), Qwtif0(MAXTOPIF,MAXQPT),
     &            Qptif1(MAXTOPIF,4,MAXQPT), Qwtif1(MAXTOPIF,MAXQPT)
        integer ::    nint(MAXTOP),           nintb(MAXTOP),
     &                nintif0(MAXTOPIF),       nintif1(MAXTOPIF),
     &                ngauss,                 ngaussb,        ngaussif,
     &                intp,
     &                maxnint
c
      end module intpt_m
