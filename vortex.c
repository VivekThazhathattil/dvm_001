#include "vortex.h"
#include <math.h>
#include <stdio.h>
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

void vortex_debug_print_position(vortex_t *v) {
  fprintf(stdout, "DEBUG: vortex: Vortex Position: (%f, %f)\n", creal(v->z), cimag(v->z));
}

void vortex_debug_print_strength(vortex_t *v) {
  fprintf(stdout, "DEBUG: vortex: Vortex Strength: %f\n", v->gamma);
}
double complex vortex_get_velocity(vortex_t *v, double complex z, double c) {
  double complex real_vort_dz, image_vort_dz, real_vort_contrib, image_vort_contrib;
  if (cabs(z - v->z) < EPSILON) {
    return 0.0 + 0.0 * I;
  }
  real_vort_dz = z - v->z;
  image_vort_dz = z - (c*c / conj(v->z));
  real_vort_contrib = (I * v->gamma) * 1 / (2.0 * M_PI * real_vort_dz);
  image_vort_contrib = -(I * v->gamma) * 1 / (2.0 * M_PI * image_vort_dz);
  return real_vort_contrib + image_vort_contrib;
}

void vortex_advect(vortex_t* v, double complex vel, double dt) {
  v->z += vel * dt;
  #if DEBUG_OLD
    vortex_debug_print_position(v);
    vortex_debug_print_strength(v);
  #endif
}

void vortex_set_nascent_vortex_strength(vortex_t *v, double Us, double dt) {
  v->gamma = 0.5 * Us * Us * dt;
}

void vortex_set_nascent_vortex_position(vortex_t *v, double o_p_mj, double theta_s) {
    double x, y;
    x = o_p_mj * cos(PI - theta_s);
    y = o_p_mj * sin(PI - theta_s);
    v->z = x + I * y;
}

vortex_t* create_nascent_vortex(vortex_t *v, double o_p_mj, double theta_s, double Us, double dt) {
    vortex_t *v_nascent;
    v_nascent = (vortex_t *) malloc(sizeof(vortex_t));
    vortex_set_nascent_vortex_position(v_nascent, o_p_mj, theta_s);
    vortex_set_nascent_vortex_strength(v_nascent, Us, dt);
    return v_nascent;
}