/* @file main.c
**
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "jnn.h"

typedef struct {
    float std_scale;
    int corrector; // corrector, window to increase total error thresh
    int seg_dist; // distance between 2 segs to be merged as one
    int window;
    int error;
    int min_seg_len; // this will affect polyA detection too, make sure it's not too high (cutting somewhere is better than nowhere)
    int chunk_size;
    int start_chunks; // num of chunks to store before processing
} jnnv3_aparam_t;

//dRNA realtime adaptor
#define JNNV3_R9_ADAPTOR { \
    .std_scale = 0.9, \
    .corrector = 1200, \
    .seg_dist = 1800, \
    .window = 300, \
    .error = 5, \
    .min_seg_len = 4000, \
    .chunk_size = 1200, \
    .start_chunks = 6, \
} \

#define JNNV3_RNA004_ADAPTOR { \
    .std_scale = 0.3, \
    .corrector = 800, \
    .seg_dist = 1200, \
    .window = 300, \
    .error = 5, \
    .min_seg_len = 1200, \
    .chunk_size = 1600, \
    .start_chunks = 6, \
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
    int corrector; // corrector, window to increase total error thresh
    int seg_dist; // distance between 2 segs to be merged as one
    int window;
    float stall_len;
    int error;
    float offset; // offset from median for polyA lo and hi error thresholds
    float range; // half height of lo and hi error thresholds
} jnnv3_pparam_t;

// dRNA realtime polyA parameters
#define JNNV3_R9_POLYA { \
    .corrector = 50, \
    .seg_dist = 200, \
    .window = 250, \
    .stall_len = 1.0, \
    .error = 30, \
    .offset = 30, \
    .range = 20, \
} \

#define JNNV3_RNA004_POLYA { \
    .corrector = 50, \
    .seg_dist = 200, \
    .window = 250, \
    .stall_len = 1.0, \
    .error = 42, \
    .offset = 44, \
    .range = 32, \
} \

typedef struct jnnv3_pstate_s {

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
void jnnv3_pcalc_param(jnnv3_pstate_t *state, jnn_pair_t adapt, jnnv3_pparam_t param, float *sig_store, int sig_size);
void jnnv3_pcore(jnnv3_pstate_t *t, jnnv3_pparam_t param, float *chunk, int current_chunk_size);