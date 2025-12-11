#ifndef SIMULATION_H
#define SIMULATION_H

#include "discrete_vortex_method.h"
#include <complex.h>
#include <stddef.h>

typedef struct SIMULATION_s{
  double U; // free-stream velocity
  double c; // cylinder radius
  double dt; // time step
  dvm_t *dvm;
} sim_t;

sim_t* sim_create(double U, double c, double dt);
void sim_kill(sim_t *s);

#endif