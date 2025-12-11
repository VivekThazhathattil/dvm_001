#include "discrete_vortex_method.h"

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
    for (size_t i = 0; i < dvm->num_vorts; ++i) {
        free(&(dvm->vs[i]));
    }
    free(dvm->vs);
    free(dvm);
}

double complex dvm_get_complex_velocity(double complex z, dvm_t *dvm) {
    double complex vort_contrib, non_vort_contrib;
    non_vort_contrib = cv_uniform_cylinder(z, dvm->U, dvm->c);
    vort_contrib = 0.0 + 0.0 * I;

    for (size_t i = 0; i < dvm->num_vorts; ++i) {
        vort_contrib += cv_vortex_with_image(&dvm->vs[i], z, dvm->c);
    }
    return non_vort_contrib + vort_contrib;
}
