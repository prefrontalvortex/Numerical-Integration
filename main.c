//
// Created by mike on 12/12/16.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
//#include <random>
#include <quadmath.h>
#include <time.h>
#include "main.h"
#include "aux.h"
#include "argparse.h"
#include "rng.h"
#include "functions.h"


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-noreturn"

#define SPECIAL_VEC
//#define KALMAN
#define MAX_ITERS 100000

void *threaded_main(void *arg);


int NUM_THREADS;
pthread_mutex_t lock_x;
precise_t REAL_VAL;


int main(int argc, char **argv) {
    precise_t (*test_fn) (precise_t); // Variable for holding the pointer to the test function
    precise_t (*integrator) (thread_data_t *); // Integrator function

    int rc, NDIM = 3;
    int seed, randVecDims, i, n, d, th;
    long iter=0, lseed, batch_size;
    precise_t val_est;
    precise_t imin, imax;
    double d_batch_size, d_imin, d_imax;
    char *function_name = (char*) emalloc(256);
    char *integrator_name = (char*) emalloc(256);

    REAL_VAL = 1;

    uint128_t count_in = 0, count_all =0;

    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    __time_t tmp_seconds = now.tv_sec;

    //// Grab parameters from the command line arguments
    arg_t *args = NULL;
    parser_populate(&args, argc, argv);
    parse_assign_i(&NUM_THREADS, "-t", args, "1");
    parse_assign_d(&d_batch_size, "-b", args, "1e7");
    parse_assign_d(&d_imin, "-a", args, "-10");
    parse_assign_d(&d_imax, "-z", args, "10");
    parse_assign_cs(&function_name, "-f", args, "identity");
    parse_assign_cs(&integrator_name, "-i", args, "riemann");



    batch_size = (long) d_batch_size;
    imin = (precise_t) d_imin;
    imax = (precise_t) d_imax;

    if ( !strcmp(function_name, "identity")) {
        fprintf(stderr, "Integrating identity(x)\n");
        test_fn = identity;
    } else if ( !strcmp(function_name, "normal")) {
        fprintf(stderr, "Integrating normal(x) \n");
        test_fn = normal;
    } else if ( !strcmp(function_name, "normaln")) {
        fprintf(stderr, "Integrating normal(x) (no coefficient)\n");
        test_fn = normal_no_coef;
    } else {
        fprintf(stderr, "Warning! No valid function name provided. Defaulting to identity function\n");
        test_fn = identity;
    }
    if ( !strcmp(integrator_name, "riem")) {
        fprintf(stderr, "Integrator: Riemann\n");
        integrator = riemann_integrator;
    } else if ( !strcmp(integrator_name, "stoc")) {
        fprintf(stderr, "Integrator: Stochastic \n");
        integrator = stochastic_integrator;
    } else if ( !strcmp(integrator_name, "nest")) {
        fprintf(stderr, "Integrator: nested sampling\n");
        die("Function not available!\n");
    } else {
        fprintf(stderr, "Warning! No valid function name provided. Defaulting to Riemann\n");
        integrator = riemann_integrator;
    }

    fprintf(stderr, "Size of uint128: %d-bit\n Precise_t: %d-bit\n", 8*sizeof(uint128_t), 8*sizeof(precise_t));
    fprintf(stderr, " Num Iters/Batch: %ld\tNum Threads: %d\n", batch_size, NUM_THREADS);

    FILE * randfile = fopen("/dev/random", "r");

    if (randfile == NULL)   die("Fatal error! Can't open file!\n");
    lseed = rand_long_from_file(randfile);
    fclose(randfile);
    seed = (int) lseed;
    srand(seed);



    SimpleAverage *avglist = new_SimpleAverage(MAX_ITERS, NUM_THREADS);
    thread_data_t payload[NUM_THREADS];
    pthread_t threads[NUM_THREADS];
    struct drand48_data rngbuf[NUM_THREADS];
    RedbearRNG_data rngbuf2[NUM_THREADS];

    if (pthread_mutex_init(&lock_x, NULL) != 0) {
        fprintf(stderr, "\nFatal: Mutex init failed\n"), exit(EXIT_FAILURE);
    }


    for (th=0; th<NUM_THREADS; th++) {
        srand48_r(lseed + th, &(rngbuf[th]));
        srand_redbear_r( (uint128_t) (lseed + th), &(rngbuf2[th]));

        payload[th].tid = th;
        payload[th].avglist = avglist;
        payload[th].batch_size = batch_size;
        payload[th].imin = imin;
        payload[th].imax = imax;
        payload[th].rngbuf_rb = &(rngbuf2[th]);
        payload[th].test_fn = test_fn;
        payload[th].integrator = integrator;
        payload[th].filename = emalloc(256);
        sprintf(payload[th].filename, "out/%ld_pi_avg_t%d.csv", (long) tmp_seconds, th);

    }

    for (th=0; th<NUM_THREADS; th++) {
        rc = pthread_create(&threads[th], NULL, threaded_main, (void *) &payload[th]);

    }
    for (th = 0; th < NUM_THREADS; th++){
        pthread_join(threads[th], NULL);
    }

    val_est = 0;
    for (th = 0; th < NUM_THREADS; th++){
        val_est += payload[th].val_est;
    }
    val_est /= (precise_t) NUM_THREADS;

    if (NUM_THREADS >= 2) {
        fprintf(stdout, "Threaded estimate of integral: %.13lf\n", (double) val_est);
    }

    return EXIT_SUCCESS;
}



void *threaded_main(void *arg) {
    thread_data_t *payload = (thread_data_t *) arg;
    RedbearRNG_data *rngbuf_rb = payload->rngbuf_rb;
    FILE *pifile;
    stopwatch_t time;
    int i, J, tid;
    double logDelta, seconds, std;
    tid = payload->tid;
    long iter = 0, batch_size = payload->batch_size;
    precise_t x, y, val_calc, val_est, delta, K_gain, alpha, var;
    val_est = 1.;
    alpha = 0.1;

//    pifile = fopen(payload->filename, "w");
//    fprintf(pifile, "thread,iters,alpha,pi_batch,pi_kalman,error,logError\n");
//    fclose(pifile);

    for (J=0; J < 1; J++) {
        startTimer(&time);

        val_calc = payload->integrator(payload);
//        val_calc = riemann_integrator(payload->test_fn, payload->imin, payload->imax, batch_size, rngbuf_rb);
//        val_est = (1-alpha)*val_est + alpha * val_calc;
        pthread_mutex_lock(&lock_x);
        val_est = simple_average_observe(payload->avglist, val_calc);
        pthread_mutex_unlock(&lock_x);

        var = payload->avglist->var;
        std = sqrt((double) var);
        delta = val_est - REAL_VAL;
        logDelta = fabs((double) delta);
        logDelta = log10(logDelta);
//        pifile = fopen(payload->filename, "a");
//        fprintf(pifile, "%d,%ld,%le,%.16lf,%.16lf,%.16lf,%.8lf\n", tid, batch_size*iter, (double) alpha,
//                (double) val_calc, (double) val_est, (double) delta, logDelta);
//        fclose(pifile);
        if (tid == 0) {
//            fprintf(stdout, "<%d>Iter: %5ld alpha: %le  sec/thr: %5.2f\n", tid, iter, (double) alpha, seconds/ (double) NUM_THREADS);
            fprintf(stdout, "Integral Calculated :  %.12lf\n", (double) val_calc);
//            fprintf(stdout, "<%d>Pi IIR  :  %.12lf\tStd: %0.12lf\n", tid, (double) val_est, std);
//            fprintf(stdout, "<%d>Diff    : %+.12lf\n", tid, (double) (val_calc - REAL_VAL));
//            fprintf(stdout, "<%d>Diff IIR: %+.12lf\t logD: %+.5lf\n", tid, (double) delta, logDelta);

        }

        iter++;
        alpha = 1./ log((double) 1+iter);
    }
    payload->val_est = val_est;

    pthread_exit(NULL);
}

//precise_t riemann_integrator(precise_t (*test_fn)(precise_t), precise_t imin, precise_t imax, long batch_size,
//                             RedbearRNG_data *rngbuf_rb) {
precise_t riemann_integrator(thread_data_t *params) {

    precise_t (*test_fn)(precise_t) = params->test_fn;
    precise_t imin  = params->imin;
    precise_t imax  = params->imax;
    long batch_size = params->batch_size;
    RedbearRNG_data *rngbuf_rb = params->rngbuf_rb;
    long iter;
    precise_t result, x, y, range, dx, accu, step;
    range = imax - imin;
    dx = range/(precise_t) batch_size;

    accu = 0;
    for (iter=0; iter<batch_size; iter++) {
        step = (precise_t) iter;
        x = (step+0.5)*dx + imin;
        y = test_fn(x);
        accu += y*dx;
    }
    return accu;
}

precise_t stochastic_integrator(thread_data_t *params) {

    precise_t (*test_fn)(precise_t) = params->test_fn;
    precise_t imin  = params->imin;
    precise_t imax  = params->imax;
    long batch_size = params->batch_size;
    RedbearRNG_data *rngbuf_rb = params->rngbuf_rb;
    long iter;
    precise_t result, x, y, range, dx, accu, step;
    range = imax - imin;
    dx = range/(precise_t) batch_size;

    accu = 0;
    for (iter=0; iter<batch_size; iter++) {
        x = rand_redbear_uniform_r(rngbuf_rb);
        x = range*x + imin;
        y = test_fn(x);
        accu += y*dx;
    }
    return accu;
}

SimpleAverage *new_SimpleAverage(long max_iters, int num_threads) {
    SimpleAverage *avglist = (SimpleAverage*) emalloc(sizeof(SimpleAverage));
//    avglist->samples = (precise_t*) emalloc(max_iters * num_threads * sizeof(precise_t));
    avglist->iters = (long*) emalloc(num_threads * sizeof(long));
    avglist->idx = 0;
    avglist->total = 0;
    avglist->vartotal = 0;
    avglist->REAL_PI = 3.1415926535897932384626433832795028841971; // for diagnostic purposes

    return avglist;
}


precise_t simple_average_observe(SimpleAverage *avglist, precise_t obs) {
    // store priors
//    avglist->mean = simple_average_mean(avglist);
    avglist->obs = obs;
    avglist->total += obs;
    precise_t delta = (obs - avglist->REAL_PI);
    avglist->vartotal += (delta * delta);
    avglist->idx++;
    avglist->error = obs - avglist->mean;
    avglist->mean = avglist->total / (precise_t) avglist->idx;
    avglist->var = avglist->vartotal  / (precise_t) avglist->idx;
    return avglist->mean;
}

#pragma clang diagnostic pop