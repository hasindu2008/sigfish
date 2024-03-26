/* @file main.c
**
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "rjnn.h"
#include "error.h"
#include "stat.h"

jnnv3_astate_t *init_jnnv3_astate(jnnv3_aparam_t param){
    jnnv3_astate_t *state = (jnnv3_astate_t *)malloc(sizeof(jnnv3_astate_t));
    MALLOC_CHK(state);

    state->prev = 0;  // previous string
    state->err = 0;       // total error
    state->prev_err = 0;  // consecutive error
    state->c = 0;         // counter
    state->w = param.corrector;       // window to increase total error thresh (corrector)
    state->start = 0;     // start pos
    state->end = 0;       // end pos
    state->adapter_found = 0;
    state->sig_length = 0;
    state->median = 0;
    state->stdev = 0;
    state->top = 0;

    state->seg_i = 0;
    state->seg_c = SIGTK_SIZE;
    state->segs = (jnn_pair_t *)malloc(sizeof(jnn_pair_t)*state->seg_c);
    MALLOC_CHK(state->segs);

    return state;
}

void free_jnnv3_astate(jnnv3_astate_t *state){
    free(state->segs);
    free(state);
    return;
}

void reset_jnnv3_astate(jnnv3_astate_t *state, jnnv3_aparam_t param){
    state->prev = 0;  // previous string
    state->err = 0;       // total error
    state->prev_err = 0;  // consecutive error
    state->c = 0;         // counter
    state->w = param.corrector;       // window to increase total error thresh (corrector)
    state->start = 0;     // start pos
    state->end = 0;       // end pos
    state->adapter_found = 0;
    state->sig_length = 0;
    state->median = 0;
    state->stdev = 0;
    state->top = 0;
    state->seg_i = 0;
}


void jnnv3_acalc_param(jnnv3_astate_t *s, jnnv3_aparam_t param, float *sig_store, int sig_size){
    s->median = medianf(sig_store ,sig_size);
    // use this with outlier rejection to fix s->stdev thresholds
    s->stdev = stdvf(sig_store,sig_size);
    s->top = s->median + (s->stdev * param.std_scale);
}

void jnnv3_acore(jnnv3_astate_t *s, jnnv3_aparam_t param, float *chunk, int current_chunk_size){

    const int seg_dist = param.seg_dist;  // distance between 2 segs to be merged as one
    const int window = param.window; //Minimum segment window size to be detected
    const int error = param.error; //Allowable error in segment algorithm
    const int min_seg_len = param.min_seg_len; //Minimum length of a segment to be constructed


    for (int i=0; i< current_chunk_size; i++){
        s->sig_length++;
        float a = chunk[i];
        //fprintf(stderr,"%d: %f\n", s->sig_length, a);
        if (a < s->top){ // If datapoint is within range
            if (!s->prev){
                s->start = s->sig_length;
                s->prev = 1;
            }
            s->c++; // increase counter
            s->w ++; // increase window corrector count
            if (s->prev_err){
                s->prev_err = 0;
            }
            if (s->c >= window && s->c >= s->w  && !(s->c % s->w )){ // if current window longer than detect limit, and corrector, and is divisible by corrector
                s->err -= 1; // drop current error count by 1
            }
        }
        else{
            if (s->prev && s->err < error){
                s->c++;
                s->err++;
                s->prev_err++;
                if (s->c >= window && s->c >= s->w  && !(s->c % s->w )){
                    s->err -= 1;
                }
            }
            else if (s->prev && s->c >= window){
                s->end = s->sig_length - s->prev_err; // go back to where error stretch began for accurate cutting
                s->prev = 0;
                if (s->seg_i && s->start - s->segs[s->seg_i-1].y < seg_dist){ // if segs very close, merge them
                    s->segs[s->seg_i-1].y = s->end;
                }
                else{
                    if(s->seg_i>=s->seg_c){
                        s->seg_c *= 2;
                        s->segs = (jnn_pair_t *)realloc(s->segs,sizeof(jnn_pair_t)*s->seg_c);
                        MALLOC_CHK(s->segs);
                    }
                    s->segs[s->seg_i].x = s->start;
                    s->segs[s->seg_i].y = s->end;
                    s->seg_i++;
                }
                s->c = 0;
                s->err = 0;
                s->prev_err = 0;
            }
            else if (s->prev) {
                s->prev = 0;
                s->c = 0;
                s->err = 0;
                s->prev_err = 0;
            }
            else if (s->seg_i && (s->segs[s->seg_i-1].y-s->segs[s->seg_i-1].x >= min_seg_len) && s->sig_length - s->segs[s->seg_i-1].y > seg_dist){
                //fprintf(stderr,"Break point: %ld; s->segs %d\n",s->sig_length, s->seg_i);
                s->prev = 0;
                s->c = 0;
                s->err = 0;
                s->prev_err = 0;
                s->adapter_found = 1;
                break;
            }
            else{
                continue;
            }
        }
    }
}

jnnv3_pstate_t *init_jnnv3_pstate(jnnv3_pparam_t param){
    jnnv3_pstate_t *state = (jnnv3_pstate_t *)malloc(sizeof(jnnv3_pstate_t));
    MALLOC_CHK(state);

    state->prev = 0;  // previous string
    state->err = 0;       // total error
    state->prev_err = 0;  // consecutive error
    state->c = 0;         // counter
    state->w = param.corrector;       // window to increase total error thresh (corrector)
    state->start = 0;     // start pos
    state->end = 0;       // end pos

    state->polya_found = 0;
    state->sig_length = 0;

    state->mean = 0;
    state->top = 0;
    state->bot = 0;

    state->seg_i = 0;
    state->seg_c = SIGTK_SIZE;
    state->segs = (jnn_pair_t *)malloc(sizeof(jnn_pair_t)*state->seg_c);
    MALLOC_CHK(state->segs);

    return state;
}

void free_jnnv3_pstate(jnnv3_pstate_t *state){
    free(state->segs);
    free(state);
    return;
}

void reset_jnnv3_pstate(jnnv3_pstate_t *state, jnnv3_pparam_t param){
    state->prev = 0;  // previous string
    state->err = 0;       // total error
    state->prev_err = 0;  // consecutive error
    state->c = 0;         // counter
    state->w = param.corrector;       // window to increase total error thresh (corrector)
    state->start = 0;     // start pos
    state->end = 0;       // end pos
    state->polya_found = 0;
    state->sig_length = 0;

    state->mean = 0;

    state->top = 0;
    state->bot = 0;
    state->seg_i = 0;
}


void jnnv3_pcalc_param(jnnv3_pstate_t *state, jnn_pair_t adapt, jnnv3_pparam_t param, float *sig_store, int sig_size) {
    jnn_pair_t p = adapt;
    ASSERT(p.y > 0);
    ASSERT(p.y < sig_size);
    state->mean = meanf(&sig_store[p.x],p.y-p.x);
    state->top = state->mean + param.offset + param.range;
    state->bot = state->mean + param.offset - param.range;
}


void jnnv3_pcore(jnnv3_pstate_t *t, jnnv3_pparam_t param, float *chunk, int current_chunk_size){


    const int seg_dist = param.seg_dist; // distance between 2 segs to be merged as one
    const int window = param.window;
    const int error = param.error;
    const float stall_len = param.stall_len;

    for(int i=0; i<current_chunk_size; i++){
        t->sig_length++;
        float a = chunk[i];
        a = a>OUTLIER_MAX ? OUTLIER_MAX : a;
        a = a<OUTLIER_MIN ? OUTLIER_MIN : a;
        if (a < t->top && a > t->bot){ // If datapoint is within range
            if (!t->prev){
                t->start = t->sig_length;
                t->prev = 1;
            }
            t->c++; // increase counter
            t->w++; // increase window corrector count
            if (t->prev_err){
                t->prev_err = 0;
            }
            if (t->c >= window && t->c >= t->w &&  !(t->c % t->w)){ // if current window longer than detect limit, and corrector, and is divisible by corrector
                t->err--; // drop current error count by 1
            }
        }
        else{
            if (t->prev && t->err < error){
                t->c++;
                t->err++;
                t->prev_err++;
                if (t->c >= window && t->c >= t->w && !(t->c % t->w)){
                    t->err--;
                }
            }
            else if (t->prev && (t->c >= window || (!t->seg_i && t->c >= window * stall_len))){
                t->end = t->sig_length - t->prev_err; // go back to where error stretch began for accurate cutting
                t->prev = 0;
                if (t->seg_i && t->start - t->segs[t->seg_i-1].y < seg_dist){ // if segs very close, merge them
                    t->segs[t->seg_i-1].y = t->end;
                }
                else{
                    if(t->seg_i>=t->seg_c){
                        t->seg_c *= 2;
                        t->segs = (jnn_pair_t *)realloc(t->segs,sizeof(jnn_pair_t)*t->seg_c);
                    }
                    t->segs[t->seg_i].x = t->start;
                    t->segs[t->seg_i].y = t->end;
                    t->seg_i++;
                }
                t->c = 0;
                t->err = 0;
                t->prev_err = 0;
            }
            else if (t->prev){
                t->prev = 0;
                t->c = 0;
                t->err = 0;
                t->prev_err = 0;
            }
        }
    }

    if(t->seg_i){
        t->polya_found = 1;
    }

}

