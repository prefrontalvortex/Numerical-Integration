//
// Created by mike on 13/12/16.
//

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#include "aux.h"

void die(const char *format, ...) {
    // Print message and Kill process
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stderr, format, argptr);
    va_end(argptr);
//    fprintf(stderr, "%s x_x", msg);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void *emalloc(size_t numBytes){ /* Error checked malloc */
    void *buffer = malloc(numBytes);
    if(!buffer){
        fprintf(stderr, "Fatal error: malloc error. Out of memory. \n  x_X\n");
        exit(EXIT_FAILURE);
    }
    return buffer;
}

void startTimer(stopwatch_t *stopwatch) {
    clock_gettime(CLOCK_REALTIME, &(stopwatch->startElapsed));
}

double getElaspedTime(stopwatch_t *stopwatch) {
    struct timespec now, startElapsed;
    startElapsed = stopwatch->startElapsed;
    clock_gettime(CLOCK_REALTIME, &now);
    __time_t tmp_seconds = (now.tv_sec - startElapsed.tv_sec);
    double seconds = (double) tmp_seconds;
    __time_t tmp_nano = (now.tv_nsec - startElapsed.tv_nsec);
    double nano = tmp_nano / 1e9;

    if (nano > 1e9) {
        tmp_nano = (startElapsed.tv_nsec - now.tv_nsec);
        nano = tmp_nano / 1e9;
    }
    return seconds + nano;
}

void checkpoint_stdprintf(int ccode, const char *format, ...) {
//    FILE *file = fopen(filename, "a");
//    if (file != NULL) {
    fprintf(stderr, "\x1B[38;5;%dm", ccode);

    va_list argptr;
    va_start(argptr, format);
    vfprintf(stderr, format, argptr);
    va_end(argptr);
    fprintf(stderr, "\x1B[0m"); // reset code

    fflush(stderr);
//    }
}