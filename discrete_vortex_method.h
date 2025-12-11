#ifndef DISCRETE_VORTEX_METHOD_H
#define DISCRETE_VORTEX_METHOD_H

#include "complex_velocity.h"
#include "vortex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define EPSILON 1e-10

typedef struct DVM_s{
  double U; // free-stream velocity
  double c; // cylinder radius
  size_t num_vorts; // number of vortices
  vortex_t *vs; // array of vortices
  double dt; // time step
} dvm_t;

dvm_t* dvm_create(double U, double c, size_t num_vorts);
void dvm_kill(dvm_t* dvm);
double complex dvm_get_complex_velocity(double complex z, dvm_t *dvm);

#endif
