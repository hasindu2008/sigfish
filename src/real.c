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
#include "jnn.h"


int16_t **get_chunks(const int16_t *raw, int64_t nsample, int chunk_size, int num_chunks){
    int16_t **chunks = (int16_t **)malloc(sizeof(int16_t *)*num_chunks);
    MALLOC_CHK(chunks);
    for (int chunk_i=0; chunk_i<num_chunks; chunk_i++){
        chunks[chunk_i] = (int16_t *)malloc(sizeof(int16_t)*chunk_size);
        MALLOC_CHK(chunks[chunk_i]);
        int cur_chunk_size = (chunk_i == num_chunks-1) ? nsample - chunk_i*chunk_size : chunk_size;
        for(int j=0; j<cur_chunk_size; j++){
            chunks[chunk_i][j] = raw[chunk_i*chunk_size + j];
            assert(chunk_i*chunk_size + j < nsample);
        }
    }
    return chunks;
}


void jnn_v3(const int16_t *raw, int64_t nsample){


    // now feed algorithm with chunks of signal simulating real-time
    int chunk_size = 1200;
    int num_chunks = (nsample + chunk_size -1) / chunk_size;

    int16_t **chunks = get_chunks(raw, nsample, chunk_size, num_chunks);

    int start_chunks = 6; //num of chunks to store before processing

    int16_t *sig_store = (int16_t *)malloc(sizeof(int16_t)*nsample);
    MALLOC_CHK(sig_store);
    int sig_store_i = 0;

    int64_t sig_length = 0;

    // this is the algo. Simple yet effective
    int8_t prev = 0;  // previous string
    int err = 0;       // total error
    int prev_err = 0;  // consecutive error
    int c = 0;         // counter
    int w = 1200;       // window to increase total error thresh
    int seg_dist = 1800;  // distance between 2 segs to be merged as one
    int start = 0;     // start pos
    int end = 0;       // end pos
    int8_t adapter_found = 0;
    int window = 300; //Minimum segment window size to be detected
    float std_scale = 0.9; //Scale factor of STDev about median
    int error = 5; //Allowable error in segment algorithm
    int min_seg_len = 4000; //Minimum length of a segment to be constructed

    int seg_c = SIGTK_SIZE;
    jnn_pair_t * segs = (jnn_pair_t *)malloc(sizeof(jnn_pair_t)*seg_c);
    MALLOC_CHK(segs);
    int seg_i = 0;

    //float mean = 0;
    float median = 0;
    float stdev = 0;
    float top = 0;
    int16_t a = 0;

    for (int chunk_i=0; chunk_i < num_chunks; chunk_i++){
        //fprintf(stderr,"processing chunk: %d\n", chunk_i);
        // print("chunk {}".format(chunk))
        int current_chunk_size = (chunk_i == num_chunks-1) ? nsample - chunk_i*chunk_size : chunk_size;

        if (chunk_i < start_chunks){
            for(int j=0; j<current_chunk_size; j++){
                sig_store[sig_store_i] = chunks[chunk_i][j];
                sig_store_i++;
                assert(sig_store_i < nsample);
            }
            continue;
        }
        fprintf(stderr,"size_i: %d\n", sig_store_i);

        if (sig_store_i > 0){
            for(int j=0; j<current_chunk_size; j++){
                sig_store[sig_store_i] = chunks[chunk_i][j];
                sig_store_i++;
                assert(sig_store_i < nsample);
            }
            int sig_size = sig_store_i;
            sig_store_i = 0;
            median = mediani16(sig_store,sig_size);
            // use this with outlier rejection to fix stdev thresholds
            stdev = stdvi16(sig_store,sig_size);
            top = median + (stdev * std_scale);
        }


        for (int i=0; i< current_chunk_size; i++){
            sig_length++;
            a = chunks[chunk_i][i];
            if (a < top){ // If datapoint is within range
                if (!prev){
                    start = sig_length;
                    prev = 1;
                }
                c++; // increase counter
                w++; // increase window corrector count
                if (prev_err){
                    prev_err = 0;
                }
                if (c >= window && c >= w && !(c % w)){ // if current window longer than detect limit, and corrector, and is divisible by corrector
                    err -= 1; // drop current error count by 1
                }
            }
            else{
                if (prev && err < error){
                    c++;
                    err++;
                    prev_err++;
                    if (c >= window && c >= w && !(c % w)){
                        err -= 1;
                    }
                }
                else if (prev && c >= window){
                    end = sig_length - prev_err; // go back to where error stretch began for accurate cutting
                    prev = 0;
                    if (seg_i && start - segs[seg_i-1].y < seg_dist){ // if segs very close, merge them
                        segs[seg_i-1].y = end;
                    }
                    else{
                        if(seg_i>=seg_c){
                            seg_c *= 2;
                            segs = (jnn_pair_t *)realloc(segs,sizeof(jnn_pair_t)*seg_c);
                            MALLOC_CHK(segs);
                        }
                        segs[seg_i].x = start;
                        segs[seg_i].y = end;
                        seg_i++;
                    }
                    c = 0;
                    err = 0;
                    prev_err = 0;
                }
                else if (prev) {
                    prev = 0;
                    c = 0;
                    err = 0;
                    prev_err = 0;
                }
                else if (seg_i && (segs[seg_i-1].y-segs[seg_i-1].x >= min_seg_len) && sig_length - segs[seg_i-1].y > seg_dist){
                    //fprintf(stderr,"Break point: %ld\n",sig_length);
                    prev = 0;
                    c = 0;
                    err = 0;
                    prev_err = 0;
                    adapter_found = 1;
                    break;
                }
                else{
                    continue;
                }
            }
        }
        if (adapter_found){
            break;
        }
    }

    if(seg_i<=0){
        assert(adapter_found == 0);
        printf(".\t.\n");
    } else{
        printf("%ld\t%ld\n",segs[0].x,segs[0].y);
    }

    for (int i=0; i<num_chunks; i++){
        free(chunks[i]);
    }
    free(chunks);
    free(segs);
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

    while((ret = slow5_get_next(&rec,sp)) >= 0){
        printf("%s\t%ld\t",rec->read_id,rec->len_raw_signal);
        jnn_v3(rec->raw_signal, rec->len_raw_signal);
    }

    if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
        fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
        exit(EXIT_FAILURE);
    }

    slow5_rec_free(rec);

    slow5_close(sp);

    return 0;
}
