#include "simulation.h"
#include <math.h>
#include <complex.h>

sim_t* sim_create(double U, double c, double dt) {
    sim_t *s;
    s = (sim_t *) malloc(sizeof(sim_t));
    s->U = U;
    s->c = c;
    s->dt = dt;
    s->dvm = dvm_create(U, c, 0);
    return s;
}

void sim_kill(sim_t *s) {
    dvm_kill(s->dvm);
    free(s);
}

int main() {
    sim_t *s;
    s = sim_create(1.0, 1.0, 0.125);
    return 0;
}