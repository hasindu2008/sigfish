
#include "sigfish.h"
#include "cdtw.h"
#include "cudtw.cuh"
#include "assert.h"

#define MALLOC_CHK(ret) { \
    if ((ret) == NULL) { \
        fprintf(stderr,"Could not allocate memory."); \
    } \
}

aln_t *init_aln2(){
    aln_t *aln = (aln_t *)malloc(sizeof(aln_t)*SECONDARY_CAP);
    MALLOC_CHK(aln);
    float score = INFINITY;
    float score2 = INFINITY;
    int32_t pos = -1;
    int32_t rid = -1;
    char d = 0;
    for (int l=0; l<SECONDARY_CAP;l++) {
        aln_t tmp = {rid,pos,pos,score,score2,d,0};
        aln[l] = tmp;
    }

    return aln;
}


void update_aln2(aln_t* aln, float score, int32_t rid, int32_t pos, char d, float *cost, int32_t qlen, int32_t rlen){
    int l=0;
    for(; l<SECONDARY_CAP; l++){
        if (score > aln[l].score){
            break;
        } else {
            continue;
        }
    }

    if(l!=0){
        for(int m=0;m<l-1;m++){
            aln[m] = aln[m+1];
        }
        aln[l-1].score = score;
        aln[l-1].pos_end = pos;
        aln[l-1].rid = rid;
        aln[l-1].d = d;

        aln[l-1].pos_st  = pos - qlen + 1;

    }
}

void dtw_single2(core_t* core,db_t* db, int32_t i) {

    if(db->slow5_rec[i]->len_raw_signal>0 && db->et[i].n>0){ //some checks to see if a good read

        aln_t *aln=init_aln2(); //initialise a alignment struct

        int64_t start_idx = db->qstart[i];  //starting index of the query
        int64_t end_idx = db->qend[i];      //ending index of the query

        int32_t qlen = end_idx - start_idx; //query chunk length

        int8_t rna = core->opt.flag & SIGFISH_RNA; // if data is RNA

        float *query = (float *)malloc(sizeof(float)*qlen);
        MALLOC_CHK(query);

        for(int j=0;j<qlen;j++){
            if (!(core->opt.flag & SIGFISH_INV) && rna){ //id rna we must reverse the events
                query[qlen-1-j] = db->et[i].event[j+start_idx].mean;
            }
            else{
                query[j] = db->et[i].event[j+start_idx].mean;
            }
        }

        for(int j=0;j<core->ref->num_ref;j++){

            int32_t rlen =core->ref->ref_lengths[j];
            float *cost = (float *)malloc(sizeof(float) * qlen * rlen);
            MALLOC_CHK(cost);

            subsequence(query, core->ref->forward[j], qlen , rlen, cost);
            for(int k=(qlen-1)*rlen; k< qlen*rlen; k+=qlen){
                float min_score = INFINITY;
                int32_t min_pos = -1;
                for(int m=0; m<qlen && k+m<qlen*rlen; m++){
                    if(cost[k+m] < min_score){
                        min_score = cost[k+m];
                        min_pos = m+k;
                    }
                }
                update_aln2(aln, min_score, j, min_pos-(qlen-1)*rlen, '+', cost, qlen, rlen);
            }


            if (!rna) { //if DNA we must consider the reverse strand as well
                subsequence(query, core->ref->reverse[j], qlen , rlen, cost);

                for(int k=(qlen-1)*rlen; k< qlen*rlen; k+=qlen){
                    float min_score = INFINITY;
                    int32_t min_pos = -1;
                    for(int m=0; m<qlen && k+m<qlen*rlen; m++){
                        if(cost[k+m] < min_score){
                            min_score = cost[k+m];
                            min_pos = m+k;
                        }
                    }
                    update_aln2(aln, min_score, j, min_pos-(qlen-1)*rlen, '-', cost, qlen, rlen);
                }

            }

            free(cost);

        }

        free(query);

        db->aln[i].score = aln[SECONDARY_CAP-1].score;
        db->aln[i].score2 = aln[SECONDARY_CAP-2].score;
        db->aln[i].pos_st = aln[SECONDARY_CAP-1].d == '+' ? aln[SECONDARY_CAP-1].pos_st : core->ref->ref_lengths[aln[SECONDARY_CAP-1].rid] - aln[SECONDARY_CAP-1].pos_end  ;
        db->aln[i].pos_end = aln[SECONDARY_CAP-1].d == '+' ? aln[SECONDARY_CAP-1].pos_end : core->ref->ref_lengths[aln[SECONDARY_CAP-1].rid] - aln[SECONDARY_CAP-1].pos_st  ;

        db->aln[i].pos_st += core->ref->ref_st_offset[aln[SECONDARY_CAP-1].rid];
        db->aln[i].pos_end += core->ref->ref_st_offset[aln[SECONDARY_CAP-1].rid];
        db->aln[i].rid = aln[SECONDARY_CAP-1].rid;
        db->aln[i].d = aln[SECONDARY_CAP-1].d;

        int mapq=(int)round(500*(db->aln[i].score2-db->aln[i].score)/db->aln[i].score);
        if(mapq>60){
            mapq=60;
        }
        db->aln[i].mapq = mapq;

        free(aln);
    }

}


void dtw_cuda_db(core_t *core, db_t *db){

    //copy from RAM to GPU

    //call GPU kernel
    int32_t i=0;
    for (i = 0; i < db->n_rec; i++) {
        dtw_single2(core,db,i);
    }
    //Copy results back

    return;
}