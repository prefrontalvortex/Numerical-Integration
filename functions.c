//
// Created by mike on 12/14/16.
//

#include <math.h>
#include "functions.h"
#include "aux.h"

#define ONE_OVER_SQRT_2PI 0.3989422804014327

precise_t identity(precise_t x) {
    return x;
}


precise_t normal(precise_t x) {
    // Standard normal, m = 0, sig=1
    return ONE_OVER_SQRT_2PI * exp( (double)(-(x*x)/2));
}

precise_t normal_no_coef(precise_t x) {
    // Standard normal, m = 0, sig=1
    return exp( (double)(-(x*x)/2));
}

