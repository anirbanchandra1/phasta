#include <cstdio>
#include <cassert>
//#include "simmetrix.h"

#ifdef __cplusplus
extern"C"{
#endif
void sim_get_pos_on_surf (double dx, double dy, double dz, int id) {
  printf("This is a dummy function. Please re-compile at phastaChef level.\n");
  printf("disp = %f,%f,%f; ID = %d\n", dx,dy,dz,id);
//  assert(0);
}

#ifdef __cplusplus
}
#endif

