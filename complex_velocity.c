#include "complex_velocity.h"
#include <math.h>
#include <complex.h>

double complex cv_uniform_cylinder(double complex z, double U, double c) {
  return -U * (1.0 - (c * c) / (z * z));
}

double complex cv_vortex_with_image(vortex_t *v, double complex z, double c) {
  return vortex_get_velocity(v, z, c);
}