/* @file main.c
**
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <slow5/slow5.h>
#include "error.h"
#include "stat.h"
#include "ref.h"
#include "sigfish.h"
#include "misc.h"
#include "jnn.h"

float **get_chunks(const float *raw, int64_t nsample, int chunk_size, int num_chunks){
    float **chunks = (float **)malloc(sizeof(float *)*num_chunks);
    MALLOC_CHK(chunks);
    for (int chunk_i=0; chunk_i<num_chunks; chunk_i++){
        chunks[chunk_i] = (float *)malloc(sizeof(float)*chunk_size);
        MALLOC_CHK(chunks[chunk_i]);
        int cur_chunk_size = (chunk_i == num_chunks-1) ? nsample - chunk_i*chunk_size : chunk_size;
        for(int j=0; j<cur_chunk_size; j++){
            chunks[chunk_i][j] = raw[chunk_i*chunk_size + j];
            assert(chunk_i*chunk_size + j < nsample);
        }
    }
    return chunks;
}

typedef struct {
    float std_scale;
    int corrector; //corrector, window to increase total error thresh
    int seg_dist; // distance between 2 segs to be merged as one
    int window;
    int error;
    int min_seg_len;
} jnnv3_param_t;

//dRNA realtime adaptor
#define JNNV3_ADAPTOR { \
    .std_scale = 0.9, \
    .corrector = 1200, \
    .seg_dist = 1800, \
    .window = 300, \
    .error = 5, \
    .min_seg_len = 4000, \
} \

typedef struct {

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

} jnnv3_state_t;


jnnv3_state_t *init_jnnv3_state(jnnv3_param_t param){
    jnnv3_state_t *state = (jnnv3_state_t *)malloc(sizeof(jnnv3_state_t));
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

void free_jnnv3_state(jnnv3_state_t *state){
    free(state->segs);
    free(state);
    return;
}

void reset_jnnv3_state(jnnv3_state_t *state, jnnv3_param_t param){
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


void jnnv3_calc_param(jnnv3_state_t *s, jnnv3_param_t param, float *sig_store, int sig_size){
    s->median = medianf(sig_store ,sig_size);
    // use this with outlier rejection to fix s->stdev thresholds
    s->stdev = stdvf(sig_store,sig_size);
    s->top = s->median + (s->stdev * param.std_scale);
}

void jnnv3_core(jnnv3_state_t *s, jnnv3_param_t param, float *chunk, int current_chunk_size){

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

void jnn_v3(const float *raw, int64_t nsample, jnnv3_param_t param, jnnv3_state_t *s){

    // now feed algorithm with chunks of signal simulating real-time
    const int chunk_size = 1200;
    const int start_chunks = 6; //num of chunks to store before processing
    const int num_chunks = (nsample + chunk_size -1) / chunk_size;
    float **chunks = get_chunks(raw, nsample, chunk_size, num_chunks);

    float *sig_store = (float *)malloc(sizeof(float)*nsample);
    MALLOC_CHK(sig_store);
    int sig_store_i = 0;


    // this is the algo. Simple yet effective
    reset_jnnv3_state(s,param);

    for (int chunk_i=0; chunk_i < num_chunks; chunk_i++){

        // print("chunk {}".format(chunk))
        int current_chunk_size = (chunk_i == num_chunks-1) ? nsample - chunk_i*chunk_size : chunk_size;
        float *chunk = chunks[chunk_i];
        //fprintf(stderr,"processing chunk: %d, size %d\n", chunk_i, current_chunk_size);

        if (chunk_i < start_chunks){ //append first few chunks to sig_store
            for(int j=0; j<current_chunk_size; j++){
                sig_store[sig_store_i] = chunk[j];
                sig_store_i++;
                assert(sig_store_i <= nsample);
            }
            //fprintf(stderr,"chunk_i: %d added to sig_store\n", chunk_i);
            continue;
        }
        //fprintf(stderr,"size_i: %d\n", sig_store_i);

        if (sig_store_i > 0){ //append one more chunk to sig_store
            for(int j=0; j<current_chunk_size; j++){
                sig_store[sig_store_i] = chunk[j];
                sig_store_i++;
                assert(sig_store_i <= nsample);
            }
            int sig_size = sig_store_i;
            sig_store_i = 0;
            jnnv3_calc_param(s, param, sig_store, sig_size);
            //fprintf(stderr,"s->median: %f, s->stdev: %f, s->top: %f\n", s->median, s->stdev, s->top);
            chunk = sig_store;
            current_chunk_size = sig_size;
        }

        jnnv3_core(s, param, chunk, current_chunk_size);

        if (s->adapter_found){
            break;
        }
    }

    if(s->seg_i<=0){
        assert(s->adapter_found == 0);
        printf(".\t.\n");
    } else{
        printf("%ld\t%ld\n",s->segs[0].x,s->segs[0].y);
    }

    for (int i=0; i<num_chunks; i++){
        free(chunks[i]);
    }
    free(chunks);
    free(sig_store);

}

int real_main(int argc, char* argv[]){

    if(argc < 2){
        fprintf(stderr,"Usage: %s <file>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    slow5_file_t *sp = slow5_open(argv[1],"r");
    if(sp==NULL){
       fprintf(stderr,"Error in opening file\n");
       exit(EXIT_FAILURE);
    }
    slow5_rec_t *rec = NULL;
    int ret=0;

    printf("read_id\tlen_raw_signal\tadapt_start\tadapt_end\n");

    const jnnv3_param_t param = JNNV3_ADAPTOR;
    jnnv3_state_t *s= init_jnnv3_state(param);

    while((ret = slow5_get_next(&rec,sp)) >= 0){
        printf("%s\t%ld\t",rec->read_id,rec->len_raw_signal);
        float *signal = signal_in_picoamps(rec);
        jnn_v3(signal, rec->len_raw_signal, param, s);
        free(signal);
    }

    if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
        fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
        exit(EXIT_FAILURE);
    }

    free_jnnv3_state(s);
    slow5_rec_free(rec);
    slow5_close(sp);

    return 0;
}
