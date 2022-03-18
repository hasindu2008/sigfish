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
    float m = mean(x, n);
    for (i = 0; i < n; i++) {
        sum += (x[i] - m) * (x[i] - m);
    }
    return sqrt(sum / n);
}

float *rolling_window(float *x, int n, int w) {
    int i = 0;

    float *t = (float*)malloc(sizeof(float)*(n-w));
    MALLOC_CHK(t);        
    for(int i=0; i<n-w; i++){
        t[i] = 0.0;
    }

    for(int i=0;i<n-w;i++){
        for(int j=0;j<w;j++){
            t[i] += x[i+j];
        }
        t[i]/=w;
    }

    return t;

}


int stats_main(int argc, char* argv[]) {

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


    if (argc-optind<1 ||  fp_help == stdout) {
        fprintf(fp_help,"Usage: sigfish stats reads.slow5 .. \n");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -h                         help\n");
        fprintf(fp_help,"   --version                  print version\n");
        //fprintf(fp_help,"   --rna                      the dataset is direct RNA\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    printf("read_id\tnsample\tmean\tstd\tw_mean\tw_std\n");

    for(int i=optind ; i<argc ;i++){

        slow5_file_t *sp = slow5_open(argv[i],"r");
        if(sp==NULL){
            ERROR("Error in opening file %s\n",argv[i]);
            exit(EXIT_FAILURE);
        }
        slow5_rec_t *rec = NULL;
        int ret=0;

        while((ret = slow5_get_next(&rec,sp)) >= 0){
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

            printf("%f\t%f\t%f\t%f",m,s,m_t,s_t);
            // for(uint64_t i=0;i<len_raw_signal;i++){
            //     double pA = TO_PICOAMPS(rec->raw_signal[i],rec->digitisation,rec->offset,rec->range);
            //     printf("%f ",pA);
            // }

            free(t);
            free(current);
            printf("\n");
        }

        if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
            fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
            exit(EXIT_FAILURE);
        }

        slow5_rec_free(rec);

        slow5_close(sp);

    }


    return 0;
}
