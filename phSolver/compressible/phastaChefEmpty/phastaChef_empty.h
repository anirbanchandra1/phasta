#ifndef PHASTA_CHEF_EMPTY_H
#define PHASTA_CHEF_EMPTY_H

#include<stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void core_get_pos_on_surf (double dx[], double dy[], double dz[], int numnp,
                           double px[], double py[], double pz[]);

void core_is_in_closure (int e_dim, int e_tag, int t_dim, int t_tag, int answer);

void core_measure_mesh (double x1[], double x2[], double x3[], int numnp,
                        double minq);

#ifdef __cplusplus
}
#endif

#endif
