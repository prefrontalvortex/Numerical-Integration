//
// Created by mike on 13/12/16.
//

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "aux.h"
#include "rng.h"

//const uint128_t REDBEAR_MAX = (uint128_t) 0xFFFFFFFFFFFFFFFFFFFFFFFF; // 2**96 -1
//const uint128_t REDBEAR_MULT = (uint128_t) 26409387504754779197847983445;
//const uint128_t  REDBEAR_INC = (uint128_t) 6094474039558795199503380795;

// knuth MMIX
const uint128_t REDBEAR_MAX = (uint128_t) 0xffffffffffffffff; // 2**64 -1
const uint128_t REDBEAR_MULT = (uint128_t) 6364136223846793005;
const uint128_t  REDBEAR_INC = (uint128_t) 1442695040888963407;

precise_t rand_uniform() {
    return rand_uniform62();
}

precise_t rand_uniform_r(struct drand48_data *rngbuffer) {
    double result;
    drand48_r(rngbuffer, &result);
    return result;
}

precise_t rand_uniform32(long rand_int) {

    return (precise_t)rand_int / (precise_t)RAND_MAX;

}


precise_t rand_from_file(FILE *randstream) {
    long rand_int;
//    fgets(&rand_int, sizeof(rand_int), randstream);
    fread(&rand_int, sizeof(rand_int), 1, randstream);

    return (precise_t)rand_int / (precise_t)RAND_MAX;

}

long rand_long_from_file( FILE * randstream) {
    long rand_int;
//    fgets(&rand_int, sizeof(rand_int), randstream);
    fread(&rand_int, sizeof(rand_int), 1, randstream);

    return rand_int;

}

precise_t rand_uniform62() {
    unsigned long accu = (unsigned long) rand();
    accu <<= 31;
    accu |= (unsigned long) rand();
    accu &= RAND62_MAX;
    return accu / (precise_t) RAND62_MAX;
}

precise_t rand_redbear_uniform_r(RedbearRNG_data *rngbuffer) {
    return (precise_t ) rand_redbear_r(rngbuffer) / (precise_t) REDBEAR_MAX;
}

uint128_t rand_redbear_r(RedbearRNG_data *rngbuffer) {
    // Combine LCG and Xor shift register
    // General LCG: Xn+1 = (aXn + c)mod m   n>= 0
    uint128_t x = rngbuffer->state_xor;
    rngbuffer->newstate = REDBEAR_MULT*rngbuffer->state + REDBEAR_INC; //
    rngbuffer->state = rngbuffer->newstate & REDBEAR_MAX; // modulo 2**96
    x ^= x >> 12; // a
    x ^= x << 25; // b
    x ^= x >> 27; // c
    rngbuffer->state_xor = x;
    x *= 0x2545F4914F6CDD1D;
    x <<= 64; // dump this into the higher bitspace
    return (rngbuffer->state ^ x) & REDBEAR_MAX;
}

uint128_t srand_redbear_r(uint128_t seed, RedbearRNG_data *rngbuffer) {
    rngbuffer->state = seed;
    rngbuffer->state_xor = seed ^ 0xDEADBEEFACE1CA75;
}
