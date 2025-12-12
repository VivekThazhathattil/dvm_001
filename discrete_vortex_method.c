#include "discrete_vortex_method.h"
#include <time.h>

dvm_t* dvm_create(double U, double c, size_t num_vorts) {
    dvm_t *dvm;
    dvm = (dvm_t *)malloc(sizeof(dvm_t));
    dvm->U = U;
    dvm->c = c;
    dvm->num_vorts = num_vorts;
    if (num_vorts > 0) {
        dvm->vs = (vortex_t *) malloc(sizeof(vortex_t) * num_vorts);
    } else {
        dvm->vs = NULL;
    }
    return dvm;
}

void dvm_kill(dvm_t* dvm) {
    //for (size_t i = 0; i < dvm->num_vorts; ++i) {
    //    fprintf(stdout, "Freeing vortex %zu\n", i);
    //    free(&(dvm->vs[i]));
    //}
    free(dvm->vs);
    free(dvm);
}

double complex dvm_get_complex_velocity(double complex z, dvm_t *dvm) {
    double complex vort_contrib, non_vort_contrib;
    non_vort_contrib = cv_uniform_cylinder(z, dvm->U, dvm->c);
    vort_contrib = 0.0 + 0.0 * I;

    for (size_t i = 0; i < dvm->num_vorts; ++i) {
        vort_contrib += cv_vortex_with_image(&(dvm->vs[i]), z, dvm->c);

        #if DEBUG_OLD
            if(isnan(creal(vort_contrib)) || isnan(cimag(vort_contrib))) {
                fprintf(stdout,\
                     "DEBUG: dvm: NaN detected in vortex contribution (vort_idx: %zu) at point (%f, %f)\n",\
                     i, creal(z), cimag(z));
            }
            if(isinf(creal(vort_contrib)) || isinf(cimag(vort_contrib))) {
                fprintf(stdout,\
                     "DEBUG: dvm: Inf detected in vortex contribution (vort_idx: %zu) at point (%f, %f)\n",\
                     i, creal(z), cimag(z));
            }
            if(isnan(creal(non_vort_contrib)) || isnan(cimag(non_vort_contrib))) {
                fprintf(stdout,\
                     "DEBUG: dvm: NaN detected in non-vortex contribution (vort_idx: %zu) at point (%f, %f)\n",\
                     i, creal(z), cimag(z));
            }
            if(isinf(creal(non_vort_contrib)) || isinf(cimag(non_vort_contrib))) {
                fprintf(stdout,\
                     "DEBUG: dvm: Inf detected in non-vortex contribution (vort_idx: %zu) at point (%f, %f)\n",\
                     i, creal(z), cimag(z));
            }
        #endif
    }

    
    return non_vort_contrib + vort_contrib;
}

void dvm_append_vortex(dvm_t* dvm, vortex_t* v) {
    dvm->vs = (vortex_t *) realloc(dvm->vs, sizeof(vortex_t) * (dvm->num_vorts + 1));
    dvm->vs[dvm->num_vorts] = *v;
    dvm->num_vorts += 1;
}

void dvm_debug_put_arbitrary_vortices(dvm_t *dvm, size_t num_vorts, \
    double x_min, double x_max, double y_min, double y_max) {
    srand((unsigned int)time(NULL));
    size_t i;
    double x, y, strength;
    vortex_t* v;
    for (i = 0; i < num_vorts; ++i) {
        v = (vortex_t *) malloc(sizeof(vortex_t));
        x = x_min + (x_max - x_min) * rand() / RAND_MAX;
        y = y_min + (y_max - y_min) * rand() / RAND_MAX;
        strength = -12.0 + 24.0 * rand() / RAND_MAX;  // strength in range [-12, 12]
        
        v->z = x + I * y;
        v->gamma = strength;
        #if DEBUG_OLD
            fprintf(stdout, "DEBUG: dvm: Debug Vortex %zu: Position = (%f, %f), Strength = %f\n", \
                i, x, y, strength);
        #endif
        dvm_append_vortex(dvm, v);
    }
}

void dvm_advect_vortices(dvm_t *dvm, double dt) {
    size_t i;
    double complex vel;
    for (i = 0; i < dvm->num_vorts; ++i) {
        vel = dvm_get_complex_velocity(dvm->vs[i].z, dvm);
        #if DEBUG_OLD
            fprintf(stdout, "DEBUG: dvm: Advecting Vortex %zu: Velocity = (%f, %f)\n", i, creal(vel), cimag(vel));
        #endif
        vortex_advect(&(dvm->vs[i]), vel, dt);
    }
}

void dvm_create_nascent_vortices() {

}