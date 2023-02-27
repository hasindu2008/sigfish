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
#include "cdtw.h"

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

} jnnv3_astate_t;


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

typedef struct {

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


void jnnv3_pcalc_param(jnnv3_pstate_t *state, jnn_pair_t adapt, float *sig_store, int sig_size){
    jnn_pair_t p = adapt;
    assert(p.y > 0);
    assert(p.y < sig_size);
    state->mean = meanf(&sig_store[p.x],p.y-p.x);
    state->top = state->mean+30+20;
    state->bot = state->mean+30-20;
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

#define QUERY_SIZE_EVENTS 250
#define QUERY_SIZE_SIG 6000
void normalise_events(event_t *rawptr,int64_t start_idx,int64_t end_idx);
aln_t *init_aln();
void update_aln(aln_t* aln, float score, int32_t rid, int32_t pos, char d, float *cost, int32_t qlen, int32_t rlen);
void update_min(int32_t *min_pos_p, float *min_score_p, float *cost, int32_t qlen, int32_t rlen, int32_t k);
void update_best_aln(aln_t *best, aln_t* aln, refsynth_t *ref);
void print_aln(int64_t start_event_idx, int64_t end_event_idx, event_table et, aln_t aln, refsynth_t *ref,  char *read_id, uint64_t len_raw_signal);
static aln_t map(refsynth_t *ref, float *raw, int64_t nsample, int polyend, char *read_id){
    assert(ref != NULL);
    assert(raw != NULL);
    assert(nsample > 0);
    assert(nsample-polyend >= QUERY_SIZE_SIG);
    int8_t rna = 1;

    aln_t best_aln = {0};
    best_aln.pos_st = -1;

    event_table et = getevents(nsample, raw, rna);
    if(et.n > 0){
        int64_t start_idx = -1;
        int64_t end_idx = -1;
        int i = 0;
        while(i < et.n && et.event[i].start < (uint64_t)polyend) i++;
        start_idx = i;
        assert((uint64_t)start_idx < et.n);
        end_idx = start_idx + QUERY_SIZE_EVENTS;

        if (start_idx + 25 > et.n ){
            fprintf(stderr,"WARNING: not enough events to map - a weird read (<25 events in %ld samples)\n",nsample-polyend);
            start_idx = 0;end_idx = 0;
        } else if(end_idx > et.n){
            fprintf(stderr,"WARNING: Only %ld events in %ld samples\n",et.n-start_idx,nsample-polyend);
            end_idx = et.n;
        }
        normalise_events(et.event,start_idx,end_idx);

        aln_t *aln=init_aln();
        int32_t qlen = end_idx - start_idx;

        float *query = (float *)malloc(sizeof(float)*qlen);
        MALLOC_CHK(query);

        for(int j=0;j<qlen;j++){
            query[qlen-1-j] = et.event[j+start_idx].mean;
        }

        for(int j=0;j<ref->num_ref;j++){
            int32_t rlen =ref->ref_lengths[j];
            float *cost = (float *)malloc(sizeof(float) * qlen * rlen);
            MALLOC_CHK(cost);
            subsequence(query, ref->forward[j], qlen , rlen, cost);
            for(int k=(qlen-1)*rlen; k< qlen*rlen; k+=qlen){
                float min_score = INFINITY;
                int32_t min_pos = -1;
                update_min(&min_pos, &min_score, cost, qlen, rlen, k);
                update_aln(aln, min_score, j, min_pos-(qlen-1)*rlen, '+', cost, qlen, rlen);
            }
            free(cost);
        }

        free(query);
        update_best_aln(&best_aln, aln, ref);
        free(aln);

        if(best_aln.pos_st > 0){
            print_aln(start_idx, end_idx, et, best_aln,  ref, read_id, nsample);
        }
    }
    free(et.event);

    return best_aln;
}

void jnn_v3(const float *raw, int64_t nsample, jnnv3_aparam_t param, jnnv3_astate_t *s, jnnv3_pparam_t pparam, jnnv3_pstate_t *t, refsynth_t *ref, char *read_id){

    // now feed algorithm with chunks of signal simulating real-time
    const int chunk_size = 1200;
    const int start_chunks = 6; //num of chunks to store before processing
    const int num_chunks = (nsample + chunk_size -1) / chunk_size;
    float **chunks = get_chunks(raw, nsample, chunk_size, num_chunks);

    float *sig_store = (float *)malloc(sizeof(float)*nsample);
    MALLOC_CHK(sig_store);
    int sig_store_i = 0;

    // this is the algo. Simple yet effective
    reset_jnnv3_astate(s,param);
    reset_jnnv3_pstate(t,pparam);
    aln_t best_aln = {0};
    best_aln.pos_st = -1;

    for (int chunk_i=0; chunk_i < num_chunks; chunk_i++){

        // print("chunk {}".format(chunk))
        int current_chunk_size = (chunk_i == num_chunks-1) ? nsample - chunk_i*chunk_size : chunk_size;
        float *chunk = chunks[chunk_i];
        //fprintf(stderr,"processing chunk: %d, size %d\n", chunk_i, current_chunk_size);

        for(int j=0; j<current_chunk_size; j++){
            sig_store[sig_store_i] = chunk[j];
            sig_store_i++;
            assert(sig_store_i <= nsample);
        }

        if (chunk_i < start_chunks){ //nothing
            continue;
            //fprintf(stderr,"chunk_i: %d added to sig_store\n", chunk_i);

        }
        //fprintf(stderr,"size_i: %d\n", sig_store_i);
        if (chunk_i == start_chunks){ //enough chunks arrived
            int sig_size = sig_store_i;
            //sig_store_i = 0;
            jnnv3_acalc_param(s, param, sig_store, sig_size);
            //fprintf(stderr,"s->median: %f, s->stdev: %f, s->top: %f\n", s->median, s->stdev, s->top);
            chunk = sig_store;
            current_chunk_size = sig_size;
        }

        if (!s->adapter_found){
            jnnv3_acore(s, param, chunk, current_chunk_size);
            if (s->adapter_found){
                jnn_pair_t p = s->segs[0];
                jnnv3_pcalc_param(t,p,sig_store,sig_store_i);
                chunk = &sig_store[p.y];
                current_chunk_size = sig_store_i-p.y;

            }
        }

        if(s->adapter_found && !t->polya_found){
            jnnv3_pcore(t, pparam,chunk,current_chunk_size);
        }

        if(t->polya_found){
            assert(s->adapter_found == 1);
            assert(t->seg_i > 0);
            assert(s->seg_i > 0);
            jnn_pair_t polya = t->segs[0];
            jnn_pair_t adapt = s->segs[0];
            int st = polya.y+adapt.y-1;
            int leftover = sig_store_i - st;
            if(leftover >= QUERY_SIZE_SIG){
                //fprintf(stderr,"leftover: %d, running DTW\n", leftover);
                if(ref){
                    best_aln=map(ref, sig_store, sig_store_i, st, read_id);
                }
                break;
            } else {
                //fprintf(stderr,"leftover: %d, waiting for more\n", leftover);
            }

        }

    }

    if(ref==NULL){
        if(s->seg_i<=0){
            assert(s->adapter_found == 0);
            printf(".\t.\t");
        } else{
            printf("%ld\t%ld\t",s->segs[0].x,s->segs[0].y);
        }
        if(t->seg_i <= 0){
            assert(t->polya_found == 0);
            printf(".\t.\n");
        } else {
            jnn_pair_t polya = t->segs[0];
            assert(polya.y + s->segs[0].y < sig_store_i);
            printf("%ld\t%ld\n",polya.x+s->segs[0].y-1, polya.y+s->segs[0].y-1);
        }
    }


    for (int i=0; i<num_chunks; i++){
        free(chunks[i]);
    }
    free(chunks);
    free(sig_store);

}



int real_main(int argc, char* argv[]){

    if(argc < 2){
        fprintf(stderr,"Usage: sigfish %s <file.blow5>\n", argv[0]);
        fprintf(stderr,"Usage: sigfish %s <file.blow5> <ref.fa>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    slow5_file_t *sp = slow5_open(argv[1],"r");
    if(sp==NULL){
       fprintf(stderr,"Error in opening file\n");
       exit(EXIT_FAILURE);
    }
    slow5_rec_t *rec = NULL;
    int ret=0;

    refsynth_t *ref = NULL;
    if(argc == 3){
        const char *ref_name = argv[2];
        assert(ref_name != NULL);
        model_t *pore_model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER); //4096 is 4^6 which is hardcoded now
        MALLOC_CHK(pore_model);
        uint32_t kmer_size = set_model(pore_model, MODEL_ID_RNA_NUCLEOTIDE);
        uint32_t flag = 0;
        flag |= SIGFISH_RNA;
        int32_t query_size = QUERY_SIZE_EVENTS;
        ref = gen_ref(ref_name, pore_model, kmer_size, flag, query_size);
        free(pore_model);

    }

    if(ref == NULL) printf("read_id\tlen_raw_signal\tadapt_start\tadapt_end\tpolya_start\tpolya_end\n");

    const jnnv3_aparam_t param = JNNV3_ADAPTOR;
    jnnv3_astate_t *s= init_jnnv3_astate(param);
    const jnnv3_pparam_t pparam = JNNV3_POLYA;
    jnnv3_pstate_t *t = init_jnnv3_pstate(pparam);

    while((ret = slow5_get_next(&rec,sp)) >= 0){
        if(ref == NULL) printf("%s\t%ld\t",rec->read_id,rec->len_raw_signal);
        float *signal = signal_in_picoamps(rec);
        jnn_v3(signal, rec->len_raw_signal, param, s, pparam, t, ref, rec->read_id);
        free(signal);
    }

    if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
        fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
        exit(EXIT_FAILURE);
    }

    free_jnnv3_astate(s);
    free_jnnv3_pstate(t);
    slow5_rec_free(rec);
    slow5_close(sp);

    if(ref != NULL){
        free_ref(ref);
    }

    return 0;
}
