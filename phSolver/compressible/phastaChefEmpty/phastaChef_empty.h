#ifndef PHASTA_CHEF_EMPTY_H
#define PHASTA_CHEF_EMPTY_H

#include<stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void sim_get_pos_on_surf (double dx[], double dy[], double dz[], int numnp,
                          double px[], double py[], double pz[]);

void sim_is_in_closure (int e_dim, int e_tag, int t_dim, int t_tag, int answer);

#ifdef __cplusplus
}
#endif

#endif
