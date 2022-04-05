/* @file  events_main.c
**
** @@
******************************************************************************/

#include "sigfish.h"
#include "misc.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


void print_events(char *rid, event_table et){

    uint32_t j = 0;
    for (j = 0; j < et.n; j++) {
        printf("%s\t%d\t%ld\t%d\t%f\t%f\n", rid, j, et.event[j].start,
                (int)et.event[j].length, et.event[j].mean,
                et.event[j].stdv);
    }
    printf("\n");
     
}

void event_func(slow5_rec_t *rec, int8_t rna){

    float *current_signal = signal_in_picoamps(rec);
    //trim(current_signal, rec->len_raw_signal);

    event_table et = getevents(rec->len_raw_signal, current_signal, rna);
    print_events(rec->read_id,et);

    free(current_signal);
    free(et.event);    
}

float mean(float *x, int n) {
    float sum = 0;
    int i = 0;
    for (i = 0; i < n; i++) {
        sum += x[i];
    }
    return sum / n;
}

float stdv(float *x, int n) {
    float sum = 0;
    int i = 0;
    float m = mean(x, n); //can reuse already calculated mean
    for (i = 0; i < n; i++) {
        sum += (x[i] - m) * (x[i] - m);
    }
    return sqrt(sum / n);
}

#define OUTLIER_MAX 1200
#define OUTLIER_MIN 0

float *rolling_window(float *x, int n, int w) {
    // int i = 0;

    float *t = (float*)malloc(sizeof(float)*(n-w));
    MALLOC_CHK(t);        
    for(int i=0; i<n-w; i++){ //can remove later with a single var that inits to 0
        t[i] = 0.0;
    }

    for(int i=0;i<n-w;i++){
        for(int j=0;j<w;j++){
            t[i] += x[i+j];
        }
        t[i]/=w;
    }

    // for(int i=0; i<n-w; i++){ //can remove later with a single var that inits to 0
    //     fprintf(stderr,"%f\t",t[i]);
    // }
    // fprintf(stderr,"\n");

    return t;
}


//adapted from
//https://github.com/Psy-Fer/deeplexicon/blob/master/scripts/dRNA_segmenter.py
float *rm_outlier(int16_t *x, int n) {
    
    float *t = (float*)malloc(sizeof(float)*(n));
    MALLOC_CHK(t); 

    for(int i=0; i<n; i++){
        if(x[i]>OUTLIER_MAX){
            t[i] =OUTLIER_MAX;
        } else if (x[i]<OUTLIER_MIN){
            t[i] = OUTLIER_MIN;
        } else {
            t[i] = x[i];
        }
    }

    return t;
}

//adapted from
//https://github.com/Psy-Fer/deeplexicon/blob/master/scripts/dRNA_segmenter.py
pair_t trim(slow5_rec_t *rec){

    int64_t nsample = rec->len_raw_signal;

    if(nsample > SIGFISH_WINDOW_SIZE){


        float *current = rm_outlier(rec->raw_signal,rec->len_raw_signal);
        float *t = rolling_window(current,nsample,SIGFISH_WINDOW_SIZE);
        free(current);
        float mn = mean(t,nsample-SIGFISH_WINDOW_SIZE);
        float std = stdv(t,nsample-SIGFISH_WINDOW_SIZE);

        //int top = (MEAN_VAL + (STDV_VAL*0.5));
        float bot = mn - (std*0.5);

        int8_t begin = 0;
        int seg_dist = 1500;
        int hi_thresh = 200000;
        int lo_thresh = 2000;
        int start=0;
        int end=0;

        pair_t segs[SIGFISH_SIZE];
        int seg_i = 0;

        int count = -1;

        for (int j=0; j<nsample-SIGFISH_WINDOW_SIZE; j++){
            float i = t[j];
            count++;
            if (i < bot && !begin){
                start = count;
                begin = 1;
            }
            else if (i < bot){
                end = count;
            }
            else if (i > bot && begin){
                if (seg_i && start - segs[seg_i-1].y < seg_dist){
                    segs[seg_i-1].y = end;
                }
                else{
                    segs[seg_i] = {start,end};    
                    seg_i++;
                }
                start = 0;
                end = 0;
                begin = 0;
            }
        }
        pair_t p = {0, 0};
        for (int i=0; i<seg_i; i++){
            int b = segs[i].y;
            int a = segs[i].x;
            if (b - a > hi_thresh){
                continue;
            }
            if (b - a < lo_thresh){
                continue;
            }
            p = {a+SIGFISH_WINDOW_SIZE/2-1, b+SIGFISH_WINDOW_SIZE/2-1};
            //fprintf(stderr,"FF %d\t%d\n",p.x, p.y);
            break;
        }

        free(t); 
        return p;   
    }
    else{
        WARNING("%s","Not enough data to trim\n");
        pair_t p = {-1,-1};
        return p;
    }

}

//adapted from https://github.com/Psy-Fer/SquiggleKit/blob/a667e461b82f0ccc0a8103a8c1759515d6e34ac9/segmenter.py
pair_t trim_polya(slow5_rec_t *rec){

    int64_t nsample = rec->len_raw_signal;

    if(nsample > 0){

        float *sig = rm_outlier(rec->raw_signal,rec->len_raw_signal);

        float mn = mean(sig,nsample);
        float std = stdv(sig,nsample);

        int top = (mn + (mn*0.5));
        float bot = mn - (std*0.75);

        int8_t prev = 0;    // previous string
        int err = 0;        // total error
        int prev_err = 0;  // consecutive error
        int c = 0;         // counter
        int w = 50;        // window to increase total error thresh
        int seg_dist = 50;  // distance between 2 segs to be merged as one
        int start = 0;     // start pos
        int end = 0 ;      // end pos
        int window = 50;
        int error = 5;
        float stall_len = 0.25;

        // segments [(start, stop)]
        pair_t segs[SIGFISH_SIZE];
        int seg_i = 0;

        for(int i=0; i<nsample; i++){
            float a = sig[i];
            if (a < top and a > bot){ // If datapoint is within range
                if (!prev){
                    start = i;
                    prev = 1;
                }
                c++; // increase counter
                w++; // increase window corrector count
                if (prev_err){
                    prev_err = 0;
                }
                if (c >= window && c >= w &&  !(c % w)){ // if current window longer than detect limit, and corrector, and is divisible by corrector
                    err--; // drop current error count by 1
                }
            }
            else{
                if (prev && err < error){
                    c++;
                    err++;
                    prev_err++;
                    if (c >= window && c >= w && !(c % w)){
                        err--;
                    }
                }
                else if (prev && (c >= window || (!seg_i && c >= window * stall_len))){
                    end = i - prev_err; // go back to where error stretch began for accurate cutting
                    prev = 0;
                    if (seg_i && start - segs[seg_i-1].y < seg_dist){ // if segs very close, merge them
                        segs[seg_i-1].y = end;
                    }
                    else{
                        segs[seg_i] = {start,end}; 
                        seg_i++;  
                    }
                    c = 0;
                    err = 0;
                    prev_err = 0;
                }
                else if (prev){
                    prev = 0;
                    c = 0;
                    err = 0;
                    prev_err = 0;
                }
            }
        }

        free(sig);

        printf("%d\t",seg_i);
        for(int i=0; i<seg_i; i++){
            if (i==seg_i-1){
                printf("%d,%d",segs[i].x,segs[i].y);
            }else{
                printf("%d,%d\t",segs[i].x,segs[i].y);
            }
        }
        printf("\n");

        //return p;   
    }
}



void stat_func(slow5_rec_t *rec, int8_t rna){
    printf("%s\t",rec->read_id);

    uint64_t len_raw_signal = rec->len_raw_signal;
    printf("%ld\t",len_raw_signal);


    float *current = signal_in_picoamps(rec);

    float m = mean(current,len_raw_signal);
    float s = stdv(current,len_raw_signal);

    float *t = rolling_window(current,len_raw_signal,SIGFISH_WINDOW_SIZE);
    int t_size = len_raw_signal-SIGFISH_WINDOW_SIZE;

    float m_t = mean(t,t_size);
    float s_t = stdv(t,t_size);
    printf("%f\t%f\t%f\t%f\t",m,s,m_t,s_t);

    pair_t p=trim(rec);
    printf("%d\t%d",p.x, p.y);
    // for(uint64_t i=0;i<len_raw_signal;i++){
    //     double pA = TO_PICOAMPS(rec->raw_signal[i],rec->digitisation,rec->offset,rec->range);
    //     printf("%f ",pA);
    // }

    free(t);
    free(current);
    printf("\n");    
}



void seg_func(slow5_rec_t *rec, int8_t rna){
    printf("%s\t",rec->read_id);
    trim_polya(rec);
}
