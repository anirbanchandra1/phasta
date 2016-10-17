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
