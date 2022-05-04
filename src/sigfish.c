/* @file sigfish.c
**
** @@
******************************************************************************/

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sigfish.h"
#include "misc.h"
#include "cdtw.h"

#include "slow5/slow5.h"
#include "../slow5lib/src/slow5_extra.h"

#include <sys/wait.h>
#include <unistd.h>


/* initialise the core data structure */
core_t* init_core(const char *fastafile, char *slow5file, opt_t opt,double realtime0) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    core->reg_list=NULL; //region list is NULL by default
    core->reg_i=0;
    core->reg_n=0;

    if(opt.region_str == NULL){
    }
    else{
        //determine if .bed
        int region_str_len = strlen(opt.region_str);
        if(region_str_len>=4 && strcmp(&(opt.region_str[region_str_len-4]),".bed")==0 ){
            STDERR("Fetching the list of regions from file: %s", opt.region_str);
            WARNING("%s", "Loading region windows from a bed file is an experimental option and not yet throughly tested.");
            WARNING("%s", "When loading windows from a bed file, output is based on reads that are unclipped. Also, there may be repeated entries when regions overlap.");
            int64_t count=0;
            char **reg_l = read_bed_regions(opt.region_str, &count);
            core->reg_list = reg_l;
            core->reg_i = 0;
            core->reg_n = count;
        }
        else{
            STDERR("Limiting to region: %s\n", opt.region_str);
        }
    }

    core->sf = slow5_open(slow5file,"r");
    if (core->sf == NULL) {
        STDERR("Error opening SLOW5 file %s\n",slow5file);
        exit(EXIT_FAILURE);
    }


    //model
    core->model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER); //4096 is 4^6 which is hardcoded now
    MALLOC_CHK(core->model);
    // core->cpgmodel = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER_METH); //15625 is 4^6 which os hardcoded now
    // MALLOC_CHK(core->cpgmodel);

    //load the model from files
    uint32_t kmer_size=0;
    //uint32_t kmer_size_meth=0;
    if (opt.model_file) {
        kmer_size=read_model(core->model, opt.model_file, MODEL_TYPE_NUCLEOTIDE);
    } else {
        if(opt.flag & SIGFISH_RNA){
            INFO("%s","builtin RNA nucleotide model loaded");
            kmer_size=set_model(core->model, MODEL_ID_RNA_NUCLEOTIDE);
        }
        else{
            kmer_size=set_model(core->model, MODEL_ID_DNA_NUCLEOTIDE);
        }
    }
    // if (opt.meth_model_file) {
    //     kmer_size_meth=read_model(core->cpgmodel, opt.meth_model_file, MODEL_TYPE_METH);
    // } else {
    //     kmer_size_meth=set_model(core->cpgmodel, MODEL_ID_DNA_CPG);
    // }
    // if( kmer_size != kmer_size_meth){
    //     ERROR("The k-mer size of the nucleotide model (%d) and the methylation model (%d) should be the same.",kmer_size,kmer_size_meth);
    //     exit(EXIT_FAILURE);
    // }
    core->kmer_size = kmer_size;


    //synthetic reference
    core->ref = gen_ref(fastafile,core->model,kmer_size,opt.flag, opt.query_size);

    core->opt = opt;

    //realtime0
    core->realtime0=realtime0;

    core->load_db_time=0;
    core->process_db_time=0;
    core->output_time=0;

    core->sum_bytes=0;
    core->total_reads=0; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)


    return core;
}

/* free the core data structure */
void free_core(core_t* core,opt_t opt) {
    free(core->model);
    // free(core->cpgmodel);

    if(core->reg_list){
        for(int64_t i=0;i<core->reg_n;i++){
            free(core->reg_list[i]);
        }
        free(core->reg_list);
    }

    free_ref(core->ref);

    slow5_close(core->sf);
    free(core);
}

/* initialise a data batch */
db_t* init_db(core_t* core) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_rec = core->opt.batch_size;
    db->n_rec = 0;

    db->mem_records = (char**)(calloc(db->capacity_rec,sizeof(char*)));
    MALLOC_CHK(db->mem_records);
    db->mem_bytes = (size_t*)(calloc(db->capacity_rec,sizeof(size_t)));
    MALLOC_CHK(db->mem_bytes);

    db->slow5_rec = (slow5_rec_t**)calloc(db->capacity_rec,sizeof(slow5_rec_t*));
    MALLOC_CHK(db->slow5_rec);

    db->current_signal = (float**)malloc(sizeof(float*) * db->capacity_rec);
    MALLOC_CHK(db->current_signal);

    db->et = (event_table*)malloc(sizeof(event_table) * db->capacity_rec);
    MALLOC_CHK(db->et);

    db->aln = (aln_t *)malloc(sizeof(aln_t) * db->capacity_rec);
    MALLOC_CHK(db->et);

    db->qstart = (int64_t *)malloc(sizeof(int64_t) * db->capacity_rec);
    MALLOC_CHK(db->qstart);

    db->qend = (int64_t *)malloc(sizeof(int64_t) * db->capacity_rec);
    MALLOC_CHK(db->qend);

    db->total_reads=0;


    return db;
}

/* load a data batch from disk */
ret_status_t load_db(core_t* core, db_t* db) {

    double load_start = realtime();

    db->n_rec = 0;
    db->sum_bytes = 0;
    db->total_reads = 0;

    ret_status_t status={0,0};
    int32_t i = 0;
    while (db->n_rec < db->capacity_rec && db->sum_bytes<core->opt.batch_size_bytes) {
        i=db->n_rec;
        db->mem_records[i] = (char *)slow5_get_next_mem(&(db->mem_bytes[i]), core->sf);

        if (db->mem_records[i] == NULL) {
            if (slow5_errno != SLOW5_ERR_EOF) {
                ERROR("Error reading from SLOW5 file %d", slow5_errno);
                exit(EXIT_FAILURE);
            }
            else {
                break;
            }
        }
        else {
            db->n_rec++;
            db->total_reads++; // candidate read
            db->sum_bytes += db->mem_bytes[i];
        }
    }

    status.num_reads=db->n_rec;
    status.num_bytes=db->sum_bytes;

    double load_end = realtime();
    core->load_db_time += (load_end-load_start);

    return status;
}


void parse_single(core_t* core,db_t* db, int32_t i){

    assert(db->mem_bytes[i]>0);
    assert(db->mem_records[i]!=NULL);
    //db->slow5_rec[i]=NULL;
    int ret=slow5_rec_depress_parse(&db->mem_records[i], &db->mem_bytes[i], NULL, &db->slow5_rec[i], core->sf);
    if(ret!=0){
        ERROR("Error parsing the record %d",i);
        exit(EXIT_FAILURE);
    }

}

void event_single(core_t* core,db_t* db, int32_t i) {

    if(db->slow5_rec[i]->len_raw_signal>0){

        int16_t* rawptr = db->slow5_rec[i]->raw_signal;
        float range = db->slow5_rec[i]->range;
        float digitisation = db->slow5_rec[i]->digitisation;
        float offset = db->slow5_rec[i]->offset;
        int32_t nsample = db->slow5_rec[i]->len_raw_signal;

        // convert to pA
        db->current_signal[i] = (float*)malloc(sizeof(float) * nsample);
        MALLOC_CHK(db->current_signal[i]);

        float raw_unit = range / digitisation;
        for (int32_t j = 0; j < nsample; j++) {
            db->current_signal[i][j] = ((float)rawptr[j] + offset) * raw_unit;
        }

        int8_t rna=0;
        if (core->opt.flag & SIGFISH_RNA){
            rna=1;
        }
        db->et[i] = getevents(nsample, db->current_signal[i], rna);


        //get the scalings
        // db->scalings[i] = estimate_scalings_using_mom(
        //     db->read[i], db->read_len[i], core->model, core->kmer_size, db->et[i]);

        // //If sequencing RNA, reverse the events to be 3'->5'
        // if (rna){
        //     event_t *events = db->et[i].event;
        //     size_t n_events = db->et[i].n;
        //     for (size_t i = 0; i < n_events/2; ++i) {
        //         event_t tmp_event = events[i];
        //         events[i]=events[n_events-1-i];
        //         events[n_events-1-i]=tmp_event;
        //     }
        // }


    }
    else{
        db->et[i].n = 0;
        db->et[i].event = NULL;
    }

}

int64_t detect_query_start(slow5_rec_t *rec, event_table et){
    int64_t start = -1;
    pair_t p=find_adaptor(rec);
    int64_t len_raw_signal = rec->len_raw_signal;
    if(p.y > 0){
        assert(p.y<len_raw_signal);

        float *current = signal_in_picoamps(rec);
        float m_a = meanf(&current[p.x],p.y-p.x);
        // float s_a = stdvf(&current[p.x],p.y-p.x);
        // float k_a = medianf(&current[p.x],p.y-p.x);

        assert(p.y > 0);
        assert(p.y < len_raw_signal);

        float *adapt_end = &current[p.y];
        pair_t polya = find_polya(adapt_end,len_raw_signal-p.y, m_a+30+20,m_a+30-20);

        uint64_t i = 0;
        if (polya.y > 0){
            assert(et.n>0);
            assert(et.event[i].start>=0);

            polya.y = polya.y + p.y;

            while(i < et.n && et.event[i].start < (uint64_t)polya.y ){
                i++;
            }
            start = i;
            if((uint64_t)start >= et.n){
                start = -1;
            }
            fprintf(stderr,"%s\tevent:%ld/%ld\traw:%ld/%ld\n",rec->read_id,start,et.n,polya.y ,len_raw_signal);
        } else {
            fprintf(stderr,"%s\t./%ld\n",rec->read_id,len_raw_signal);
        }

    }
    return start;

}

void normalise_single(core_t* core,db_t* db, int32_t i) {


    if(db->slow5_rec[i]->len_raw_signal>0 && db->et[i].n>0){

        int64_t start_idx;
        int64_t end_idx;
        int64_t n = db->et[i].n;

        int8_t from_sig_end= core->opt.flag & SIGFISH_END;

        if(!from_sig_end){ //map query start

            start_idx =  core->opt.prefix_size;
            if(core->opt.prefix_size < 0){
                start_idx = detect_query_start(db->slow5_rec[i], db->et[i]);
                if(start_idx < 0){
                    WARNING("Autodetect query start failed for %s",db->slow5_rec[i]->read_id);
                    start_idx = 0;
                    end_idx = 0;
                    db->et[i].n = 0;
                }
            }
            end_idx = start_idx+core->opt.query_size;

            if (start_idx > n) {
                WARNING("Read %s is ignored (%ld events < %d prefix)",db->slow5_rec[i]->read_id, db->et[i].n, core->opt.prefix_size);
                start_idx = 0;
                end_idx = 0;
                db->et[i].n = 0;
            }
            if(end_idx > n){
                WARNING("Read %s is too short (%ld events < %d prefix+query_size)",db->slow5_rec[i]->read_id, db->et[i].n, core->opt.prefix_size+core->opt.query_size);
                end_idx = db->et[i].n;
            }

        }
        else{ //map query end
            start_idx = db->et[i].n - core->opt.prefix_size - core->opt.query_size;
            end_idx = db->et[i].n - core->opt.prefix_size;
            if (start_idx < 0) {
                WARNING("Read %s is too short (%ld events < %d prefix+query_size)",db->slow5_rec[i]->read_id, db->et[i].n, core->opt.prefix_size+core->opt.query_size);
                start_idx = 0;
            }
            if(end_idx <0 ){
                WARNING("Read %s is ignored (%ld events < %d prefix)",db->slow5_rec[i]->read_id, db->et[i].n, core->opt.prefix_size);
                end_idx = 0;
                db->et[i].n = 0;
            }
        }
        db->qstart[i] = start_idx;
        db->qend[i] = end_idx;


        float event_mean = 0;
        float event_var = 0;
        float event_stdv = 0;
        float num_samples = end_idx-start_idx;

        event_t *rawptr = db->et[i].event;

        for(int64_t j=start_idx; j<end_idx; j++){
            event_mean += rawptr[j].mean;
        }
        event_mean /= num_samples;
        for(int64_t j=start_idx; j<end_idx; j++){
            event_var += (rawptr[j].mean-event_mean)*(rawptr[j].mean-event_mean);
        }
        event_var /= num_samples;
        event_stdv = sqrt(event_var);

        for(int64_t j=start_idx; j<end_idx; j++){
            rawptr[j].mean = (rawptr[j].mean-event_mean)/event_stdv;
        }
    }

}

aln_t *init_aln(){
    aln_t *aln = (aln_t *)malloc(sizeof(aln_t)*SECONDARY_CAP);
    MALLOC_CHK(aln);
    float score = INFINITY;
    float score2 = INFINITY;
    int32_t pos = -1;
    int32_t rid = -1;
    char d = 0;
    for (int l=0; l<SECONDARY_CAP;l++) aln[l] = {rid,pos,pos,score,score2,d,0};

    return aln;
}

void update_aln(aln_t* aln, float score, int32_t rid, int32_t pos, char d, float *cost, int32_t qlen, int32_t rlen){
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

        Path p;
        if(subsequence_path(cost, qlen, rlen, pos, &p)){
            if(p.k<=0){
                fprintf(stderr,"Could not find path as size is 0\n");
                aln[l-1].pos_st = -1;
            }else{
                aln[l-1].pos_st = p.py[0];
                if(p.py[p.k-1] != pos){
                    fprintf(stderr,"SOme shit happened in the backtracking\n");
                    fprintf(stderr,"%d %d %d %d %d\n",qlen,rlen,pos,aln[l-1].pos_st,p.py[p.k-1]);
                }

            }
            free(p.px);
            free(p.py);
        }
        else {
            fprintf(stderr,"Could not find path. %d %d %d\n",qlen,rlen,pos);
            aln[l-1].pos_st = -1;
        }

    }

}


void dtw_single(core_t* core,db_t* db, int32_t i) {

    if(db->slow5_rec[i]->len_raw_signal>0 && db->et[i].n>0){

        aln_t *aln=init_aln();

        int64_t start_idx = db->qstart[i];
        int64_t end_idx = db->qend[i];
        int64_t n =  db->et[i].n;

        int8_t from_sig_end= core->opt.flag & SIGFISH_END;
        int32_t qlen;

        if(!from_sig_end){ //map query start
            // start_idx =  core->opt.prefix_size;
            // end_idx = start_idx+core->opt.query_size;
            //qlen = end_idx > n ? n -start_idx : core->opt.query_size;
            qlen = end_idx - start_idx;
        }
        else{  //map query end
            // start_idx = n - core->opt.prefix_size - core->opt.query_size;
            // end_idx = n - core->opt.prefix_size;
            //qlen = start_idx < 0 ? end_idx : core->opt.query_size;
            qlen = end_idx - start_idx;
            assert(qlen>=0);
        }

        int8_t rna = core->opt.flag & SIGFISH_RNA;

        float *query = (float *)malloc(sizeof(float)*qlen);
        MALLOC_CHK(query);

        for(int j=0;j<qlen;j++){
            if (!(core->opt.flag & SIGFISH_INV) && rna){
                query[qlen-1-j] = db->et[i].event[j+start_idx].mean;
            }
            else{
                query[j] = db->et[i].event[j+start_idx].mean;
            }
        }

        //fprintf(stderr,"numref %d\n",core->ref->num_ref)    ;
        for(int j=0;j<core->ref->num_ref;j++){

            int32_t rlen =core->ref->ref_lengths[j];
            float *cost = (float *)malloc(sizeof(float) * qlen * rlen);
            MALLOC_CHK(cost);

            //fprintf(stderr,"%d,%d\n",qlen,rlen);

            if(!(core->opt.flag & SIGFISH_DTW)){
                // fprintf(stderr,"query: ");
                // for(int k=0;k<qlen;k++){
                //     fprintf(stderr,"%f,",query[k]);
                // }
                //                 fprintf(stderr,"\n");
                // fprintf(stderr,"Ref: ");

                // for(int k=0;k<rlen;k++){
                //     fprintf(stderr,"%f,",core->ref->forward[j][k]);
                // }
                // fprintf(stderr,"\n\n");
                subsequence(query, core->ref->forward[j], qlen , rlen, cost);
                for(int k=(qlen-1)*rlen; k< qlen*rlen; k+=qlen){
                    float min_score = INFINITY;
                    int32_t min_pos = -1;
                    for(int m=0;m<qlen && k+m<qlen*rlen;m++){
                        if(cost[k+m] < min_score){
                            min_score = cost[k+m];
                            min_pos = m+k;
                        }
                    }
                    update_aln(aln, min_score, j, min_pos-(qlen-1)*rlen, '+', cost, qlen, rlen);
                }

                // for(int k=(qlen-1)*rlen; k< qlen*rlen; k++){
                //     update_aln(aln, cost[k], j, k-(qlen-1)*rlen, '+',);
                //     // if(cost[k]<score){
                //     //     score2=score;
                //     //     score = cost[k];
                //     //     pos = k-(qlen-1)*rlen;
                //     //     rid = j;
                //     //     d = '+';
                //     // }
                // }
            }
            else{
                std_dtw(query, core->ref->forward[j], qlen , rlen, cost, 0);
                int k=qlen*rlen-1;
                update_aln(aln, cost[k], j, k-(qlen-1)*rlen, '+', cost, qlen, rlen);
                // if(cost[k]<score){
                //     score2=score;
                //     score = cost[k];
                //     pos = k-(qlen-1)*rlen;
                //     rid = j;
                //     d = '+';
                // }
            }
            // for(int k=0;k<qlen;k++){
            //     for(int l=0;l<rlen;l++){
            //         fprintf(stderr,"%f,",cost[k*rlen+l]);
            //     }
            //     fprintf(stderr,"\n");
            // }
            // fprintf(stderr,"\n");
            // exit(0);

            if (!rna) {
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
                    update_aln(aln, min_score, j, min_pos-(qlen-1)*rlen, '-', cost, qlen, rlen);
                }

                // for(int k=(qlen-1)*rlen; k< qlen*rlen; k++){
                //     update_aln(aln, cost[k], j, k-(qlen-1)*rlen, '-');
                //     // if(cost[k]<score){
                //     //     score2=score;
                //     //     score = cost[k];
                //     //     pos = k-(qlen-1)*rlen;
                //     //     rid = j;
                //     //     d = '-';
                //     // }
                // }
            }

            free(cost);

        }

        free(query);


        db->aln[i].score = aln[SECONDARY_CAP-1].score;
        db->aln[i].score2 = aln[SECONDARY_CAP-2].score;
        db->aln[i].pos_st = aln[SECONDARY_CAP-1].d == '+' ? aln[SECONDARY_CAP-1].pos_st : core->ref->ref_lengths[aln[SECONDARY_CAP-1].rid] - aln[SECONDARY_CAP-1].pos_end  ;
        db->aln[i].pos_end = aln[SECONDARY_CAP-1].d == '+' ? aln[SECONDARY_CAP-1].pos_end : core->ref->ref_lengths[aln[SECONDARY_CAP-1].rid] - aln[SECONDARY_CAP-1].pos_st  ;
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


void work_per_single_read(core_t* core,db_t* db, int32_t i){
    parse_single(core,db,i);
    event_single(core,db,i);
    normalise_single(core,db,i);
    dtw_single(core,db,i);

}

void process_db(core_t* core,db_t* db){
    double proc_start = realtime();
    work_db(core, db, work_per_single_read);
    double proc_end = realtime();
    core->process_db_time += (proc_end-proc_start);
}

/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db) {

    double output_start = realtime();

    int32_t i = 0;
    for (i = 0; i < db->n_rec; i++) {
        // printf(">%s\tLN:%d\tEVENTSTART:%d\tEVENTEND:%d\n",
        //        db->slow5_rec[i]->read_id, (int)db->et[i].n,
        //        (int)db->et[i].start, (int)db->et[i].end);
        // uint32_t j = 0;
        // for (j = 0; j < db->et[i].n; j++) {
        //     printf("{%d,%f,%f,%f}\t", (int)db->et[i].event[j].start,
        //            db->et[i].event[j].length, db->et[i].event[j].mean,
        //            db->et[i].event[j].stdv);
        // }
        // printf("\n");

        if(db->slow5_rec[i]->len_raw_signal>0 && db->et[i].n>0){
            // Output of results
            uint64_t start_idx =  db->qstart[i];
            uint64_t end_idx =  db->qend[i];
            uint64_t query_size =  end_idx-start_idx;
            float block_len = db->aln[i].pos_end - db->aln[i].pos_st;
            float residue = block_len - db->aln[i].score*block_len/query_size ;


            printf("%s\t",db->slow5_rec[i]->read_id); // Query sequence name
            printf("%d\t%ld\t%ld\t", db->et[i].n, start_idx,end_idx); // Query sequence length, start, end (in terms of events)
            printf("%c\t",db->aln[i].d); // Direction
            printf("%s\t",core->ref->ref_names[db->aln[i].rid]); // Target sequence name
            printf("%d\t",core->ref->ref_seq_lengths[db->aln[i].rid]); // Target sequence length


            printf("%d\t",db->aln[i].pos_st); // Target start
            printf("%d\t",db->aln[i].pos_end); // Target end
            printf("%d\t",(int)round(residue)); // Number of residues //todo check this
            printf("%d\t",(int)round(block_len)); //  Alignment block length //todo check this
            printf("%d\t",db->aln[i].mapq); // Mapq
            printf("tp:A:P\t");
            printf("d1:f:%.2f\t",db->aln[i].score); // distance of the best match
            printf("d2:f:%.2f\n",db->aln[i].score2); // distance of the second best matcj

        }

    }

    core->sum_bytes += db->sum_bytes;
    core->total_reads += db->total_reads;


    //core->read_index = core->read_index + db->n_rec;
    double output_end = realtime();
    core->output_time += (output_end-output_start);

}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_rec; ++i) {
        free(db->current_signal[i]);
        free(db->et[i].event);
        free(db->mem_records[i]);
    }
}

/* completely free a data batch */
void free_db(db_t* db) {

    int32_t i = 0;
    for (i = 0; i < db->capacity_rec; ++i) {
        slow5_rec_free(db->slow5_rec[i]);
    }
    free(db->slow5_rec);
    free(db->mem_records);
    free(db->mem_bytes);
    free(db->current_signal);
    free(db->qstart);
    free(db->qend);

    free(db->et);
    free(db->aln);
    //free(db->scalings);

    free(db);
}

/* initialise user specified options */
void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->batch_size = 512;
    opt->batch_size_bytes = 20*1000*1000;
    opt->num_thread = 8;
    opt->region_str = NULL; //whole genome processing if null

    opt->model_file = NULL;
    opt->meth_model_file = NULL;

    opt->debug_break=-1;

    opt->prefix_size = 50;
    opt->query_size = 250;

}
