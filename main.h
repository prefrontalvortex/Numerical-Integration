//
// Created by mike on 12/14/16.
//

#ifndef INTEGRATOR_MAIN_H
#define INTEGRATOR_MAIN_H

#include "aux.h"
#include "rng.h"



typedef struct _SimpleAverage {
    precise_t obs;          // observations (normal about x, sigma=0.1);
    precise_t *samples;     // array of samples
    precise_t total;
    precise_t vartotal;
    precise_t mean;
    precise_t var;
    precise_t error;
    precise_t REAL_PI;
    long idx;
    long *iters;             // iterations run for
} SimpleAverage;

typedef struct _thread_data_t {
    int tid;
    long batch_size;
    precise_t val_calc;
    precise_t val_est;
    precise_t (*test_fn) (precise_t);
    precise_t imin;
    precise_t imax;
    SimpleAverage *avglist;
    RedbearRNG_data *rngbuf_rb;

    char *filename;

} thread_data_t;

typedef struct _integrator_data_t {
    int tid;
    long batch_size;
    precise_t val_calc;
    precise_t val_est;
    precise_t (*test_fn) (precise_t);
    precise_t imin;
    precise_t imax;
    RedbearRNG_data *rngbuf_rb;
    char *filename;

} integrator_data_t;


SimpleAverage *new_SimpleAverage(long max_iters, int num_threads);
precise_t simple_average_observe(SimpleAverage *avglist, precise_t obs);
precise_t simple_average_mean(SimpleAverage *avglist);
//precise_t riemann_integrator(precise_t (*test_fn)(precise_t), precise_t imin, precise_t imax, long batch_size,
//                             RedbearRNG_data *rngbuf_rb);
precise_t riemann_integrator(thread_data_t *params);

#endif //INTEGRATOR_MAIN_H
