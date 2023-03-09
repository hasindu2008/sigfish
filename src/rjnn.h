/* @file main.c
**
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "jnn.h"

typedef struct {
    float std_scale;
    int corrector; //corrector, window to increase total error thresh
    int seg_dist; // distance between 2 segs to be merged as one
    int window;
    int error;
    int min_seg_len;
} jnnv3_aparam_t;

//dRNA realtime adaptor
#define JNNV3_ADAPTOR { \
    .std_scale = 0.9, \
    .corrector = 1200, \
    .seg_dist = 1800, \
    .window = 300, \
    .error = 5, \
    .min_seg_len = 4000, \
} \

typedef struct jnnv3_astate_s {

    int8_t prev;  // previous string
    int err;       // total error
    int prev_err;  // consecutive error
    int c;         // counter
    int w;       // window to increase total error thresh (corrector)
    int start;     // start pos
    int end;       // end pos
    int8_t adapter_found;
    int64_t sig_length;
    float median;
    float stdev;
    float top;

    int seg_c;
    jnn_pair_t * segs;
    int seg_i;

} jnnv3_astate_t;


jnnv3_astate_t *init_jnnv3_astate(jnnv3_aparam_t param);

void free_jnnv3_astate(jnnv3_astate_t *state);

void reset_jnnv3_astate(jnnv3_astate_t *state, jnnv3_aparam_t param);

void jnnv3_acalc_param(jnnv3_astate_t *s, jnnv3_aparam_t param, float *sig_store, int sig_size);

void jnnv3_acore(jnnv3_astate_t *s, jnnv3_aparam_t param, float *chunk, int current_chunk_size);

typedef struct {
    int corrector; //corrector, window to increase total error thresh
    int seg_dist; // distance between 2 segs to be merged as one
    int window;
    float stall_len;
    int error;
} jnnv3_pparam_t;

//dRNA realtime polyA parameters
#define JNNV3_POLYA { \
    .corrector = 50, \
    .seg_dist = 200, \
    .window = 250, \
    .stall_len = 1.0, \
    .error = 30, \
} \

typedef struct jnnv3_pstate_s{

    int8_t prev;  // previous string
    int err;       // total error
    int prev_err;  // consecutive error
    int c;         // counter
    int w;       // window to increase total error thresh (corrector)
    int start;     // start pos
    int end;       // end pos
    int8_t polya_found;
    int64_t sig_length;

    float mean;
    float top;
    float bot;

    int seg_c;
    jnn_pair_t * segs;
    int seg_i;

} jnnv3_pstate_t;


jnnv3_pstate_t *init_jnnv3_pstate(jnnv3_pparam_t param);
void free_jnnv3_pstate(jnnv3_pstate_t *state);
void reset_jnnv3_pstate(jnnv3_pstate_t *state, jnnv3_pparam_t param);
void jnnv3_pcalc_param(jnnv3_pstate_t *state, jnn_pair_t adapt, float *sig_store, int sig_size);
void jnnv3_pcore(jnnv3_pstate_t *t, jnnv3_pparam_t param, float *chunk, int current_chunk_size);