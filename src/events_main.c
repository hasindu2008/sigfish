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



static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {"version", no_argument, 0, 'V'},              //2
    {"output",required_argument, 0, 'o'},          //3 output to a file [stdout]
    {"rna",no_argument,0,0},                       //4 if RNA
    {0, 0, 0, 0}};




void print_events(char *rid, event_table et){

    printf("read_id\tevent_idx\tevent_start\tevent_len\tevent_mean\tevent_std\n");
    uint32_t j = 0;
    for (j = 0; j < et.n; j++) {
        printf("%s\t%d\t%ld\t%d\t%f\t%f\n", rid, j, et.event[j].start,
                (int)et.event[j].length, et.event[j].mean,
                et.event[j].stdv);
    }
    printf("\n");
     
}

int events_main(int argc, char* argv[]) {

    const char* optstring = "o:hV";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    int8_t rna = 0;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='V'){
            fprintf(stdout,"sigfish %s\n",SIGFISH_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c == 0 && longindex == 4){ //if RNA
            rna = 1;
        } 
    }


    if (argc-optind<2 ||  fp_help == stdout) {
        fprintf(fp_help,"Usage: sigfish events reads.slow5 read_id1 .. \n");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -h                         help\n");
        fprintf(fp_help,"   --version                  print version\n");
        fprintf(fp_help,"   --rna                      the dataset is direct RNA\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    slow5_file_t *slow5file = slow5_open(argv[optind], "r");
    if (!slow5file) {
        ERROR("cannot open %s. \n", argv[optind]);
        exit(EXIT_FAILURE);
    }

    int ret_idx = slow5_idx_load(slow5file);
    if (ret_idx < 0) {
        ERROR("Error loading index file for %s\n", argv[optind]);
        exit(EXIT_FAILURE);
    }

    slow5_rec_t *rec = NULL;

    for(int i=optind+1 ; i<argc ;i++){
        char *read_id = argv[i];
        fprintf(stderr,"Read ID %s\n",read_id);
        int ret = slow5_get(read_id, &rec, slow5file);
        if(ret < 0){
            ERROR("%s","Error when fetching the read\n");
            exit(EXIT_FAILURE);
        }

        float *current_signal = signal_in_picoamps(rec);

        trim(current_signal, rec->len_raw_signal);

        event_table et = getevents(rec->len_raw_signal, current_signal, rna);
        print_events(read_id,et);
        
        free(current_signal);
        free(et.event);

    }

    slow5_rec_free(rec);   
    slow5_idx_unload(slow5file);
    slow5_close(slow5file);

    return 0;
}
