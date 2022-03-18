/* @file  misc.c
**
** @@
******************************************************************************/

#include "sigfish.h"
#include "misc.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


pair_t trim(float *current, int32_t nsample){

    fprintf(stderr,"nsample %d\n",nsample);
    if(nsample > SIGFISH_WINDOW_SIZE){
        //int top = (MEAN_VAL + (STDV_VAL*0.5));
        int bot = (SIGFISH_MEAN_VAL - (SIGFISH_STDV_VAL*0.5));

        int8_t begin = 0;
        int seg_dist = 1500;
        int hi_thresh = 200000;
        int lo_thresh = 2000;
        int start=-1,end=1;

        pair_t segs[SIGFISH_SIZE];
        int seg_i = 0;

        int count = -1;

        float *t = (float*)malloc(sizeof(float) * (nsample-SIGFISH_WINDOW_SIZE));
        MALLOC_CHK(t);        
        for(int i=0; i<nsample-SIGFISH_WINDOW_SIZE; i++){
            t[i] = 0;
        }

        for(int i=0;i<nsample-SIGFISH_WINDOW_SIZE;i++){
            for(int j=0;j<SIGFISH_WINDOW_SIZE;j++){
                t[i] += current[i+j];
            }
            t[i]/=SIGFISH_WINDOW_SIZE;
        }
        // for(int i=0; i<nsample-SIGFISH_WINDOW_SIZE; i++){
        //     fprintf(stderr,"%f\t",t[i]);
        // }
        // fprintf(stderr,"\n");


        for (int i=0; i<nsample-SIGFISH_WINDOW_SIZE; i++){
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
        fprintf(stderr,"%d\t%d\n",count, seg_i);
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
            p = {a, b};
            fprintf(stderr,"%d\t%d\n",p.x, p.y);
            break;
        }

        free(t); 
        return p;   
    }
    else{
        fprintf(stderr,"Not enough data to trim\n");
        pair_t p = {-1,-1};
        return p;
    }

}


