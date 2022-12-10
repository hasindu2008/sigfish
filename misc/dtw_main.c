#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "cdtw.h"

void normalise(float *data, int size){
    int n = size;
    float mean = 0;
    float var = 0;
    float stdv = 0;

    for(int j= 0; j<size; j++){
        mean += data[j];
    }
    mean /= n;

    for(int j= 0; j<size; j++){
        float a = data[j]-mean;
        var += a*a;
    }

    var /= n;
    stdv =  sqrt(var);


    for(int j= 0; j<size; j++){
        data[j] = (data[j]-mean)/stdv;
    }

    return;
}

int main(int argc, char **argv){

    double t1,t2;
    t1 = clock();

    if(argc != 5){
        fprintf(stderr,"Usage: %s ref.txt query.txt ref_size query_size\n", argv[0]);
        exit(1);
    }

    FILE *ref = fopen(argv[1], "r");
    if(ref == NULL){
        fprintf(stderr, "Error: Cannot open file %s\n", argv[1]);
        exit(1);
    }
    FILE *query = fopen(argv[2], "r");
    if(query == NULL){
        fprintf(stderr, "Error: Cannot open file %s\n", argv[2]);
        exit(1);
    }

    int ref_size = atoi(argv[3]);
    int query_size = atoi(argv[4]);

    float *ref_data = malloc(sizeof(float) * ref_size);
    float *query_data = malloc(sizeof(float) * ref_size);

    for(int i = 0; i < ref_size; i++){
        assert(fscanf(ref, "%f", &ref_data[i])==1);
    }
    for(int i = 0; i < query_size; i++){
        assert(fscanf(query, "%f", &query_data[i])==1);
    }

    fclose(ref);
    fclose(query);

    normalise(ref_data, ref_size);
    normalise(query_data, query_size);

    float *cost = (float *)malloc(sizeof(float) * query_size * ref_size);

    subsequence(query_data, ref_data, query_size , ref_size, cost);

    float min_score = INFINITY;
    int min_pos = -1;
    for(int k=(query_size-1)*ref_size; k< query_size*ref_size; k+=query_size){
        for(int m=0;m<query_size && k+m<query_size*ref_size;m++){
            if(cost[k+m] < min_score){
                min_score = cost[k+m];
                min_pos = m+k;
            }
        }
    }

    fprintf(stderr, "min_score: %f, min_pos: %d\n", min_score, min_pos-(query_size-1)*ref_size);

    free(cost);
    free(ref_data);
    free(query_data);

    t2 = clock();
    fprintf(stderr, "Time: %f\n", (t2-t1)/CLOCKS_PER_SEC);

    return 0;
}