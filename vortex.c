#include "vortex.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832
#endif

vortex_t vortex_create(double complex z, double gamma) {
  vortex_t v;
  v.z = z;
  v.gamma = gamma;
  return v;
}

void vortex_set_position(vortex_t *v, double complex z) {
  v->z = z;
}

void vortex_set_gamma(vortex_t *v, double gamma) {
  v->gamma = gamma;
}

double complex vortex_get_position(vortex_t *v) {
  return v->z;
}

double vortex_get_gamma(vortex_t *v) {
  return v->gamma;
}

double complex vortex_get_velocity(vortex_t *v, double complex z, double c) {
  double complex real_vort_dz, image_vort_dz, real_vort_contrib, image_vort_contrib;
  real_vort_dz = z - v->z;
  image_vort_dz = z - (c*c / conj(v->z));
  real_vort_contrib = (I * v->gamma) * 1 / (2.0 * M_PI * real_vort_dz);
  image_vort_contrib = -(I * v->gamma) * 1 / (2.0 * M_PI * image_vort_dz);
  return real_vort_contrib + image_vort_contrib;
}
