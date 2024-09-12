/* @file sigfish.c
**
** @@
******************************************************************************/
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sigfish.h>
#include "misc.h"
#include "cdtw.h"
#include "stat.h"
#include "jnn.h"
#include "rjnn.h"
#include "str.h"

#include "slow5/slow5.h"
#include "../slow5lib/src/slow5_extra.h"

#include <sys/wait.h>
#include <unistd.h>

enum sigfish_log_level_opt _log_level = LOG_VERB;

int8_t drna_detect(slow5_file_t *sp){

    const slow5_hdr_t* hdr = sp->header;
    int8_t rna = 0;
    char *exp =slow5_hdr_get("experiment_type", 0, hdr);
    if(exp==NULL){
        WARNING("%s","experiment_type not found in SLOW5 header. Assuming genomic_dna");
        return 0;
    }
    if (strcmp(exp,"genomic_dna")==0){
        rna = 0;
    }else if (strcmp(exp,"rna")==0){
        rna = 1;
    } else {
        WARNING("Unknown experiment type: %s. Assuming genomic_dna", exp);
    }

    for(uint32_t  i=1; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("experiment_type", i, hdr);
        if (strcmp(curr, exp)){
            WARNING("Experiment type mismatch: %s != %s in read group %d. Defaulted to %s", curr, exp, i, exp);
        }
    }
    return rna;
}

int8_t pore_detect(slow5_file_t *sp) {

    const slow5_hdr_t* hdr = sp->header;
    int8_t pore = 0;
    char *kit =slow5_hdr_get("sequencing_kit", 0, hdr);
    if(kit==NULL){
        WARNING("%s","sequencing_kit not found in SLOW5 header. Assuming R9.4.1");
        return 0;
    }
    if (strstr(kit,"114")!=NULL){
        pore = OPT_PORE_R10;
    } else if (strstr(kit,"rna004")!=NULL){
        pore = OPT_PORE_RNA004;
    } else {
        pore = OPT_PORE_R9;
    }

    for(uint32_t  i=1; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("sequencing_kit", i, hdr);
        if (strcmp(curr, kit)){
            WARNING("sequencing_kit type mismatch: %s != %s in read group %d. Defaulted to %s", curr, kit, i, kit);
        }
    }
    return pore;
}


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
            VERBOSE("Fetching the list of regions from file: %s", opt.region_str);
            WARNING("%s", "Loading region windows from a bed file is an experimental option and not yet throughly tested.");
            WARNING("%s", "When loading windows from a bed file, output is based on reads that are unclipped. Also, there may be repeated entries when regions overlap.");
            int64_t count=0;
            char **reg_l = read_bed_regions(opt.region_str, &count);
            core->reg_list = reg_l;
            core->reg_i = 0;
            core->reg_n = count;
        }
        else{
            VERBOSE("Limiting to region: %s\n", opt.region_str);
        }
    }

    core->sf = slow5_open(slow5file,"r");
    if (core->sf == NULL) {
        VERBOSE("Error opening SLOW5 file %s\n",slow5file);
        exit(EXIT_FAILURE);
    }

    //drna_mismatch(core->sf, opt.flag & SIGFISH_RNA);

    if(drna_detect(core->sf)) {
        opt.flag |= SIGFISH_RNA;
        VERBOSE("%s","Detected RNA data. --rna was set automatically.");
    }

    if(opt.pore==NULL){
        int8_t pore = pore_detect(core->sf);
        if(pore){
            opt.flag |= SIGFISH_R10;
            opt.pore_flag = pore;
            if (pore == OPT_PORE_R10) VERBOSE("%s","Detected R10 data. --pore r10 was set automatically.");
            if (pore == OPT_PORE_RNA004) VERBOSE("%s","Detected RNA004 data. --pore rna004 was set automatically.");
            if (pore == OPT_PORE_R10 && (opt.flag & SIGFISH_RNA)){
                ERROR("%s","R10 RNA data does not exist! But header header indicates that the data is R10 RNA.");
                exit(EXIT_FAILURE);
            }
        }
    }

    //model
    core->model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER); //4096 is 4^6 which is hardcoded now
    MALLOC_CHK(core->model);

    //load the model from files
    uint32_t kmer_size=0;
    if (opt.model_file) {
        kmer_size=read_model(core->model, opt.model_file, MODEL_TYPE_NUCLEOTIDE);
    } else {
        if(opt.flag & SIGFISH_RNA){
            if(opt.flag & SIGFISH_R10){
                INFO("%s","builtin RNA004 nucleotide model loaded");
                kmer_size=set_model(core->model, MODEL_ID_RNA_RNA004_NUCLEOTIDE);
            } else {
                INFO("%s","builtin RNA R9 nucleotide model loaded");
                kmer_size=set_model(core->model, MODEL_ID_RNA_R9_NUCLEOTIDE);
            }
        }
        else{
            if(opt.flag & SIGFISH_R10){
                INFO("%s","builtin DNA R10 nucleotide model loaded");
                kmer_size=set_model(core->model, MODEL_ID_DNA_R10_NUCLEOTIDE);
            } else{
                INFO("%s","builtin DNA R9 nucleotide model loaded");
                kmer_size=set_model(core->model, MODEL_ID_DNA_R9_NUCLEOTIDE);
            }
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
    core->parse_time=0;
    core->event_time=0;
    core->normalise_time=0;
    core->dtw_time=0;

    core->sum_bytes=0;
    core->total_reads=0; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)

    core->prefix_fail=0;
    core->ignored=0;
    core->too_short=0;

#ifdef HAVE_ACC
    if (core->opt.flag & SIGFISH_ACC) {
        VERBOSE("%s","Initialising accelator");
    }
#endif

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

#ifdef HAVE_ACC
    if (core->opt.flag & SIGFISH_ACC) {
        VERBOSE("%s","Freeing accelator");
    }
#endif

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

    db->out = (char **)calloc(db->capacity_rec, sizeof(char *));
    MALLOC_CHK(db->out);

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
    db->prefix_fail = 0;
    db->ignored=0;
    db->too_short=0;

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

    ASSERT(db->mem_bytes[i]>0);
    ASSERT(db->mem_records[i]!=NULL);
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

int64_t detect_query_start(slow5_rec_t *rec, event_table et, int8_t pore){
    int64_t start = -1;
    jnn_pair_t p=find_adaptor(rec, pore);
    int64_t len_raw_signal = rec->len_raw_signal;
    if(p.y > 0){
        ASSERT(p.y<len_raw_signal);

        float *current = signal_in_picoamps(rec);
        float m_a = meanf(&current[p.x],p.y-p.x);
        // float s_a = stdvf(&current[p.x],p.y-p.x);
        // float k_a = medianf(&current[p.x],p.y-p.x);

        ASSERT(p.y > 0);
        ASSERT(p.y < len_raw_signal);

        float *adapt_end = &current[p.y];
        jnn_pair_t polya = find_polya(adapt_end,len_raw_signal-p.y, m_a+30+20,m_a+30-20, pore);

        uint64_t i = 0;
        if (polya.y > 0){
            ASSERT(et.n>0);
            ASSERT(et.event[i].start>=0);

            polya.y = polya.y + p.y;

            while(i < et.n && et.event[i].start < (uint64_t)polya.y ){
                i++;
            }
            start = i;
            if((uint64_t)start >= et.n){
                start = -1;
            }
            //fprintf(stderr,"%s\tevent:%ld/%ld\traw:%ld/%ld\n",rec->read_id,start,et.n,polya.y ,len_raw_signal);
        } else {
            //fprintf(stderr,"%s\t./%ld\n",rec->read_id,len_raw_signal);
        }

        free(current);

    }
    return start;

}

void normalise_events(event_t *rawptr,int64_t start_idx,int64_t end_idx){
    float event_mean = 0;
    float event_var = 0;
    float event_stdv = 0;
    float num_samples = end_idx-start_idx;

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

void normalise_single(core_t* core,db_t* db, int32_t i) {


    if(db->slow5_rec[i]->len_raw_signal>0 && db->et[i].n>0){

        int64_t start_idx;
        int64_t end_idx;
        int64_t n = db->et[i].n;

        int8_t from_sig_end= core->opt.flag & SIGFISH_END;

        if(!from_sig_end){ //map query start

            start_idx =  core->opt.prefix_size;
            if(core->opt.prefix_size < 0){
                start_idx = detect_query_start(db->slow5_rec[i], db->et[i], core->opt.pore_flag);
                if(start_idx < 0){
                    db->prefix_fail++;
                    LOG_TRACE("Autodetect query start failed for %s",db->slow5_rec[i]->read_id);
                    start_idx = 50; //fall back to 50 events start
                    // end_idx = 0;
                    // db->et[i].n = 0;
                }
            }
            end_idx = start_idx+core->opt.query_size;

            if (start_idx +25 > n ) { //min query size 25
                LOG_TRACE("Read %s is ignored (%ld events < %ld prefix)",db->slow5_rec[i]->read_id, db->et[i].n, start_idx);
                start_idx = 0;
                end_idx = 0;
                db->et[i].n = 0;
                db->ignored++;
            }
            else if(end_idx > n){
                LOG_TRACE("Read %s is too short (%ld events < %ld prefix+query_size)",db->slow5_rec[i]->read_id, db->et[i].n, start_idx+core->opt.query_size);
                end_idx = db->et[i].n;
                db->too_short++;
            }

        }
        else{ //map query end
            start_idx = db->et[i].n - core->opt.prefix_size - core->opt.query_size;
            end_idx = db->et[i].n - core->opt.prefix_size;
            if (start_idx < 0) {
                LOG_TRACE("Read %s is too short (%ld events < %d prefix+query_size)",db->slow5_rec[i]->read_id, db->et[i].n, core->opt.prefix_size+core->opt.query_size);
                start_idx = 0;
                db->too_short++;
            }
            if(end_idx <0 ){
                LOG_TRACE("Read %s is ignored (%ld events < %d prefix)",db->slow5_rec[i]->read_id, db->et[i].n, core->opt.prefix_size);
                end_idx = 0;
                db->et[i].n = 0;
                db->ignored++;
            }
        }
        db->qstart[i] = start_idx;
        db->qend[i] = end_idx;

        normalise_events(db->et[i].event,start_idx,end_idx);

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
    for (int l=0; l<SECONDARY_CAP;l++) {
        aln_t tmp = {rid,pos,pos,score,score2,d,0,NULL,0};
        aln[l] = tmp;
    }
    return aln;
}

void free_aln(aln_t *aln){
    for (int l=0; l<SECONDARY_CAP;l++) {
        free(aln[l].r2qevent_map);
    }
    free(aln);
}


static index_pair_t *path_to_map(Path p, int32_t len){

    index_pair_t *r2qevent_map = (index_pair_t *)malloc(sizeof(index_pair_t)*len);
    MALLOC_CHK(r2qevent_map);
    for(int i=0; i<len; i++){
        r2qevent_map[i].start = -1;
        r2qevent_map[i].stop = -1;
    }

    ASSERT(p.k>0);
    int ref_st = p.py[0];
    int prev_query_idx = -1;

    for(int i=0;i<p.k;i++){

        // fprintf(stderr,"%d %d, ",p.py[i],p.px[i]);

        int ref_idx = p.py[i]-ref_st;
        int query_idx = p.px[i];

        ASSERT(ref_idx<len);

        if(r2qevent_map[ref_idx].start == -1){
            r2qevent_map[ref_idx].start = query_idx;
        }
        r2qevent_map[ref_idx].stop = query_idx;

        if(prev_query_idx == query_idx){
            r2qevent_map[ref_idx].start = r2qevent_map[ref_idx].stop = -1;
        }
        prev_query_idx = query_idx;

    }
    // fprintf(stderr,"\n");
    // for(int i=0; i<len; i++){
    //     fprintf(stderr, "%d %d, ",r2qevent_map[i].start,r2qevent_map[i].stop);

    // }
    // fprintf(stderr,"\n\n");

    return r2qevent_map;
}



void update_aln(aln_t* aln, float score, int32_t rid, int32_t pos, char d, float *cost, int32_t qlen, int32_t rlen, int8_t backtrace){
    int l=0;
    for(; l<SECONDARY_CAP; l++){
        if (score > aln[l].score){
            break;
        } else {
            continue;
        }
    }

    if(l!=0){
        if(aln[0].r2qevent_size){
            aln[0].r2qevent_size=0;
            free(aln[0].r2qevent_map);
        }
        for(int m=0;m<l-1;m++){
            aln[m] = aln[m+1];
        }
        aln[l-1].score = score;
        aln[l-1].pos_end = pos;
        aln[l-1].rid = rid;
        aln[l-1].d = d;

        if(backtrace){
            Path p;
            if(subsequence_path(cost, qlen, rlen, pos, &p)){
                if(p.k<=0){
                    fprintf(stderr,"Could not find path as size is 0\n");
                    aln[l-1].pos_st = -1;
                }else{
                    aln[l-1].pos_st = p.py[0];
                    if(p.py[p.k-1] != pos){
                        fprintf(stderr,"Some shit happened in the backtracking\n");
                        fprintf(stderr,"%d %d %d %d %d\n",qlen,rlen,pos,aln[l-1].pos_st,p.py[p.k-1]);
                    }

                    int len = aln[l-1].pos_end - aln[l-1].pos_st + 1;
                    ASSERT(len >= 0);
                    aln[l-1].r2qevent_size = len;
                    aln[l-1].r2qevent_map = path_to_map(p, len);

                }

                free(p.px);
                free(p.py);
            }
            else {
                fprintf(stderr,"Could not find path. %d %d %d\n",qlen,rlen,pos);
                aln[l-1].pos_st = -1;
            }
        } else {
            aln[l-1].r2qevent_size = 0;
            aln[l-1].r2qevent_map = NULL;
            //estimate
            ASSERT(pos>=0);
            aln[l-1].pos_st = pos - qlen/2;
            aln[l-1].pos_st = aln[l-1].pos_st < 0 ? 0 : aln[l-1].pos_st;

        }

    }
}

static char *paf_str(aln_t *aln, char *read_id, char *rname, uint64_t start_raw_idx , uint64_t end_raw_idx, uint64_t query_size, int8_t rna, uint64_t len_raw_signal, uint64_t rlength){

    kstring_t str;
    kstring_t *sp = &str;
    str_init(sp, 4000);

    float block_len = aln->pos_end - aln->pos_st;
    float residue = block_len -aln->score*block_len/(query_size) ;

    // if(db->aln[i].score>70){
    //     continue;
    // }

    sprintf_append(sp, "%s\t",read_id); // read id name
    sprintf_append(sp, "%ld\t%ld\t%ld\t", len_raw_signal, start_raw_idx, end_raw_idx); // Raw signal length, start and end
    sprintf_append(sp, "%c\t",aln->d); // Direction
    sprintf_append(sp, "%s\t",rname); // reference sequence name
    sprintf_append(sp, "%d\t",rlength); // reference sequence length


    sprintf_append(sp, "%d\t",aln->pos_st); // Reference start
    sprintf_append(sp, "%d\t",aln->pos_end); // Reference end
    sprintf_append(sp, "%d\t",(int)round(residue)); // Number of residues //todo check this
    sprintf_append(sp, "%d\t",(int)round(block_len)); //  Alignment block length //todo check this
    sprintf_append(sp, "%d\t",aln->mapq); // Mapq
    sprintf_append(sp, "tp:A:P\t");
    sprintf_append(sp, "d1:f:%.2f\t",aln->score); // distance of the best match
    sprintf_append(sp, "d2:f:%.2f\t",aln->score2); // distance of the second best matcj
    sprintf_append(sp,"en:i:%d",query_size); // number of events in the mapped segment

    sprintf_append(sp, "\n");
    str.s[str.l] = '\0';
    return sp->s;
}

//todo: so inefficient as it loops through multiple times for RNA - but just to get the functionality
static char *r2qevent_map_to_ss(aln_t *aln, int64_t qstart, event_table et, int8_t rna){
    index_pair_t *base_to_event_map = aln->r2qevent_map;
    int32_t n_kmers = aln->r2qevent_size;

    if(rna){
        int end = base_to_event_map[n_kmers-1].stop;
        ASSERT(end != -1);
        // for(int i=0; i<n_kmers; i++){
        //     fprintf(stderr,"%d %d, ",base_to_event_map[i].start,base_to_event_map[i].stop);
        // }
        // fprintf(stderr,"\n");

        for(int i=0; i<n_kmers; i++){
            if(base_to_event_map[i].start != -1){
                ASSERT(base_to_event_map[i].stop != -1);
                base_to_event_map[i].start = end - base_to_event_map[i].start;
                base_to_event_map[i].stop = end - base_to_event_map[i].stop;
            }
        }

        // for(int i=0; i<n_kmers; i++){
        //     fprintf(stderr,"%d %d, ",base_to_event_map[i].start,base_to_event_map[i].stop);
        // }
        // fprintf(stderr,"\n");

    }

    for(int i=0; i<n_kmers; i++){
        if(base_to_event_map[i].start != -1){
            ASSERT(base_to_event_map[i].stop != -1);
            base_to_event_map[i].start += qstart;
            base_to_event_map[i].stop += qstart;
        }
    }

    //lazily copy pasted from https://github.com/hasindu2008/f5c/blob/master/src/resquiggle.c for now
    int64_t signal_start_point = -1;
    int64_t signal_end_point = -1;

    kstring_t str;
    kstring_t *sp = &str;
    str_init(sp, 4000);
    int64_t ci = 0; //current index
    int64_t mi = 0;
    int64_t d = 0; //deletion count
    int8_t ff = 1; //first start flag
    int matches = 0;

    if (rna){
        for (int j = 0; j < n_kmers/2; ++j) {
            index_pair_t tmp= base_to_event_map[j];
            base_to_event_map[j]=base_to_event_map[n_kmers-1-j];
            base_to_event_map[n_kmers-1-j]=tmp;
        }
        for (int j = 0; j < n_kmers; ++j) {
            int32_t tmp = base_to_event_map[j].start;
            base_to_event_map[j].start = base_to_event_map[j].stop;
            base_to_event_map[j].stop = tmp;
        }
    }

    for (int j=0;j<n_kmers; j++){
        int32_t start_event_idx = base_to_event_map[j].start;
        int32_t end_event_idx = base_to_event_map[j].stop;
        if(start_event_idx == -1){ //deletion from ref
            ASSERT(end_event_idx == -1);
            signal_start_point = signal_end_point = -1;
            if(!ff){
                ASSERT(j!=0);
                d++;
            }

        } else {
            ASSERT(end_event_idx != -1);
            //ASSERT(start_event_idx <= end_event_idx);

            signal_start_point = et.event[start_event_idx].start; //inclusive
            if(ff) {
                ff=0;
            }
            signal_end_point = et.event[end_event_idx].start + (int)et.event[end_event_idx].length; //non-inclusive

            if(d>0){
                sprintf_append(sp,"%dD",d);
                d=0;
            }
            if(j==0) ci = signal_start_point;
            ci += (mi =  signal_start_point - ci);
            ASSERT(mi>=0); //todo remove assert for performance
            if(mi) sprintf_append(sp,"%dI",(int)mi);
            ci += (mi = signal_end_point-signal_start_point);

            ASSERT(mi>=0); //todo remove assert for performance
            if(mi) {
                matches++;
                sprintf_append(sp,"%d,",(int)mi);
            }


        }
    }

    str.s[str.l] = '\0';
    return sp->s;

}

static char *sam_str(aln_t *aln, char *read_id, char *rname, uint64_t start_raw_idx , uint64_t end_raw_idx, uint64_t qlen, int64_t qstart, event_table et, int8_t rna) {
    kstring_t str;
    kstring_t *sp = &str;
    str_init(sp, 4000);

    int flag = aln->d == '+' ? 0 : 16;
    sprintf_append(sp, "%s\t%d\t", read_id, flag); //qname, flag
    sprintf_append(sp, "%s\t%ld\t%d\t", rname, (long)aln->pos_st+1, aln->mapq); //rname, pos, mapq
    sprintf_append(sp, "%ldM\t%c\t%d\t%d\t",qlen, '*', 0, 0); //cigar, rnext, pnext, tlen

    sprintf_append(sp, "%c\t%c\t",'*','*'); //seq, qual

    uint64_t post_st = rna ? aln->pos_end : aln->pos_st;
    uint64_t post_end = rna ?aln->pos_st : aln->pos_end;

    sprintf_append(sp, "si:Z:%ld,%ld,%ld,%ld\t",start_raw_idx, end_raw_idx, post_st, post_end);

    char *ss = r2qevent_map_to_ss(aln, qstart, et, rna);
    sprintf_append(sp, "ss:Z:%s",ss);
    free(ss);

    sprintf_append(sp, "\n");
    str.s[str.l] = '\0';
    return sp->s;
}

void update_min(int32_t *min_pos_p, float *min_score_p, float *cost, int32_t qlen, int32_t rlen, int32_t k){
    float min_score = INFINITY;
    int32_t min_pos = -1;
    for(int m=0;m<qlen && k+m<qlen*rlen;m++){
        if(cost[k+m] < min_score){
            min_score = cost[k+m];
            min_pos = m+k;
        }
    }
    *min_pos_p = min_pos;
    *min_score_p = min_score;
}

void update_best_aln(aln_t *best, aln_t* aln, refsynth_t *ref){

    best->score = aln[SECONDARY_CAP-1].score;
    best->score2 = aln[SECONDARY_CAP-2].score;
    best->pos_st = aln[SECONDARY_CAP-1].d == '+' ? aln[SECONDARY_CAP-1].pos_st : ref->ref_lengths[aln[SECONDARY_CAP-1].rid] - aln[SECONDARY_CAP-1].pos_end  ;
    best->pos_end = aln[SECONDARY_CAP-1].d == '+' ? aln[SECONDARY_CAP-1].pos_end : ref->ref_lengths[aln[SECONDARY_CAP-1].rid] - aln[SECONDARY_CAP-1].pos_st  ;

    best->pos_st += ref->ref_st_offset[aln[SECONDARY_CAP-1].rid];
    best->pos_end += ref->ref_st_offset[aln[SECONDARY_CAP-1].rid];
    best->rid = aln[SECONDARY_CAP-1].rid;
    best->d = aln[SECONDARY_CAP-1].d;

    int mapq=(int)round(500*(best->score2-best->score)/best->score);
    if(mapq>60){
        mapq=60;
    }
    best->mapq = mapq;
    best->r2qevent_map = aln[SECONDARY_CAP-1].r2qevent_map;
    best->r2qevent_size = aln[SECONDARY_CAP-1].r2qevent_size;

}


static void aln_to_str(core_t* core,db_t* db, int32_t i){

    if(db->slow5_rec[i]->len_raw_signal>0 && db->et[i].n>0){
        // Output of results
        uint64_t start_event_idx =  db->qstart[i];
        uint64_t end_event_idx =  db->qend[i]-1;
        ASSERT(start_event_idx>=0 && start_event_idx<db->et[i].n);
        ASSERT(end_event_idx>=0 && end_event_idx<db->et[i].n);
        uint64_t start_raw_idx = db->et[i].event[start_event_idx].start; //inclusive
        uint64_t end_raw_idx = db->et[i].event[end_event_idx].start + db->et[i].event[end_event_idx].length; //exclusive

        uint64_t query_size =  end_event_idx-start_event_idx;
        uint64_t len_raw_signal = db->slow5_rec[i]->len_raw_signal;
        uint64_t rlength = core->ref->ref_seq_lengths[db->aln[i].rid];
        char *read_id = db->slow5_rec[i]->read_id;
        char *rname = core->ref->ref_names[db->aln[i].rid];
        int8_t rna = core->opt.flag & SIGFISH_RNA;

        ASSERT(end_raw_idx <= len_raw_signal);

        if(core->opt.flag & SIGFISH_SAM){ //can bring duplicate stuff in output_db for paf here
            db->out[i] = sam_str(&db->aln[i], read_id, rname, start_raw_idx , end_raw_idx, query_size, start_event_idx,db->et[i], rna);
        } else {
            db->out[i] = paf_str(&db->aln[i], read_id, rname, start_raw_idx , end_raw_idx, query_size, rna, len_raw_signal, rlength);
        }
    } else {
        db->out[i] = NULL;
    }

}

void dtw_single(core_t* core,db_t* db, int32_t i) {

    if(db->slow5_rec[i]->len_raw_signal>0 && db->et[i].n>0){

        aln_t *aln=init_aln();

        int64_t start_idx = db->qstart[i];
        int64_t end_idx = db->qend[i];
        //int64_t n =  db->et[i].n;

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
            ASSERT(qlen>=0);
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
                    update_min(&min_pos, &min_score, cost, qlen, rlen, k);
                    update_aln(aln, min_score, j, min_pos-(qlen-1)*rlen, '+', cost, qlen, rlen, 1);
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
                update_aln(aln, cost[k], j, k-(qlen-1)*rlen, '+', cost, qlen, rlen, 1);
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
                    update_aln(aln, min_score, j, min_pos-(qlen-1)*rlen, '-', cost, qlen, rlen, 1);
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

        update_best_aln(&(db->aln[i]), aln, core->ref);

        aln_to_str(core,db,i);
        free_aln(aln);

    }

}


void work_per_single_read(core_t* core,db_t* db, int32_t i){
    parse_single(core,db,i);
    event_single(core,db,i);
    normalise_single(core,db,i);
    dtw_single(core,db,i);

}

void align_db(core_t* core, db_t* db) {
#ifdef HAVE_ACC
    if (core->opt.flag & SIGFISH_ACC) {
        VERBOSE("%s","Aligning reads with accel");
        work_db(core,db,dtw_single);
    }
#endif

    if (!(core->opt.flag & SIGFISH_ACC)) {
        //fprintf(stderr, "cpu\n");
        work_db(core,db,dtw_single);
    }
}


void process_db(core_t* core,db_t* db){
    double proc_start = realtime();

    if(core->opt.flag & SIGFISH_PRF || core->opt.flag & SIGFISH_ACC){
        double a = realtime();
        work_db(core,db,parse_single);
        double b = realtime();
        core->parse_time += (b-a);

        a = realtime();
        work_db(core,db,event_single);
        b = realtime();
        core->event_time += (b-a);

        a = realtime();
        work_db(core,db,normalise_single);
        b = realtime();
        core->normalise_time += (b-a);

        a = realtime();
        align_db(core,db);
        b = realtime();
        core->dtw_time += (b-a);
    } else {
        work_db(core, db, work_per_single_read);
    }

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
           fputs(db->out[i],stdout);
        }

    }
    fflush(stdout);

    core->sum_bytes += db->sum_bytes;
    core->total_reads += db->total_reads;
    core->prefix_fail += db->prefix_fail;
    core->ignored += db->ignored;
    core->too_short += db->too_short;


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
        free(db->out[i]); db->out[i]=NULL;
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
    free(db->out);
    //free(db->scalings);

    free(db);
}

/* initialise user specified options */
void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->batch_size = 512;
    opt->batch_size_bytes = 20*1000*1000;
    opt->pore = NULL;
    opt->pore_flag = 0;
    opt->num_thread = 8;
    opt->region_str = NULL; //whole genome processing if null

    opt->model_file = NULL;
    opt->meth_model_file = NULL;

    opt->debug_break=-1;

    opt->prefix_size = 50;
    opt->query_size = 250;

#ifdef HAVE_ACC
    opt->flag |= SIGFISH_ACC;
#endif


}


enum sigfish_log_level_opt get_log_level(){
    return _log_level;
}

void set_log_level(enum sigfish_log_level_opt level){
    _log_level = level;
}


//realtime stuff

//todo init_opt()

sigfish_state_t *init_sigfish(const char *ref_name, int num_channels, sigfish_opt_t opt){
    sigfish_state_t *state = (sigfish_state_t *)malloc(sizeof(sigfish_state_t));
    MALLOC_CHK(state);
    state->num_channels = num_channels;
    ASSERT(opt.num_thread>0 && opt.num_thread<1000);
    state->num_thread = opt.num_thread;
    state->opt = opt;
    ASSERT(opt.dtw_cutoff > 0 && opt.dtw_cutoff < 10000);
    ASSERT(opt.query_size_events > 0 && opt.query_size_events < 10000);
    ASSERT(opt.query_size_sig > 0 && opt.query_size_events < 100000);
    int size_diff = (opt.query_size_sig - opt.samples_per_event * opt.query_size_events);
    size_diff = size_diff < 0 ? -size_diff : size_diff;
    if(size_diff > 500){
        WARNING("The opt.query_size_sig - opt.samples_per_event * opt.query_size_events is a bit too large (%d). Are you sure about what you are doing?", size_diff);
    }
    float dtw_dist_diff = (opt.dtw_cutoff - 70.0/250*opt.query_size_events); //70 score for 250 seem to hold true for rna002 and rna004
    dtw_dist_diff = dtw_dist_diff < 0 ? -dtw_dist_diff : dtw_dist_diff;
    if(dtw_dist_diff > 10){
        WARNING("The opt.dtw_cutoff - 70/250*opt.query_size_events is a bit too large (%f). Are you sure about what you are doing?", dtw_dist_diff);
    }
    //do some checks

    state->status = (enum sigfish_status*)calloc(num_channels,sizeof(enum sigfish_status));
    MALLOC_CHK(state->status);
    state->reads = (sigfish_rstate_t *)calloc(num_channels,sizeof(sigfish_rstate_t));
    MALLOC_CHK(state->reads);
    state->s = (jnnv3_astate_t **)calloc(num_channels,sizeof(jnnv3_astate_t *));
    MALLOC_CHK(state->s);
    state->t = (jnnv3_pstate_t **)calloc(num_channels,sizeof(jnnv3_pstate_t *));
    MALLOC_CHK(state->t);

    jnnv3_aparam_t param = JNNV3_R9_ADAPTOR;
    jnnv3_pparam_t pparam = JNNV3_R9_POLYA;
    if (opt.pore == OPT_PORE_RNA004) {
        jnnv3_aparam_t atmp = JNNV3_RNA004_ADAPTOR;
        param = atmp;
        jnnv3_pparam_t ptmp = JNNV3_RNA004_POLYA;
        pparam = ptmp;
    }

    for(int i=0;i<num_channels;i++){
        state->reads[i].c_raw_signal = 10000;
        state->reads[i].raw_signal = (float *)malloc(sizeof(float)*10000);
        state->reads[i].read_number=-1;
        MALLOC_CHK(state->reads[i].raw_signal);
        state->s[i] = init_jnnv3_astate(param);
        state->t[i] = init_jnnv3_pstate(pparam);
        state->reads[i].read_id = NULL;

    }

    state->ref = NULL;
    if(ref_name){

        model_t *pore_model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER); //4096 is 4^6 which is hardcoded now
        MALLOC_CHK(pore_model);
        uint32_t kmer_size = set_model(pore_model, MODEL_ID_RNA_R9_NUCLEOTIDE);
        uint32_t flag = 0;
        flag |= SIGFISH_RNA;
        if(opt.no_full_ref == 0) {
            flag |= SIGFISH_REF;
        }
        int32_t query_size = state->opt.query_size_events;
        state->ref= gen_ref(ref_name, pore_model, kmer_size, flag, query_size);
        free(pore_model);

    }

    state->debug_paf = NULL;
    state->debug = NULL;
    if(opt.debug_paf){
        if ( strcmp(opt.debug_paf,"-") == 0 ){
            state->debug_paf = stdout;
        } else {
            state->debug_paf = fopen(opt.debug_paf,"w");
            F_CHK(state->debug_paf,opt.debug_paf);
        }


        state->debug = (char **)malloc(sizeof(char *)*num_channels);
        MALLOC_CHK(state->debug);
    }

    return state;
}

void free_sigfish(sigfish_state_t *state){
    for(int i=0;i<state->num_channels;i++){
        free(state->reads[i].raw_signal);
        free(state->reads[i].read_id);
        free_jnnv3_astate(state->s[i]);
        free_jnnv3_pstate(state->t[i]);
    }
    free(state->s);
    free(state->t);
    free(state->status);
    free(state->reads);
    if(state->ref) free_ref(state->ref);

    if(state->debug_paf){
        if(state->debug_paf!=stdout) fclose(state->debug_paf);
        free(state->debug);
    }

    free(state);

}

char *sprintf_aln(int64_t start_event_idx, int64_t end_event_idx, event_table et, aln_t aln, refsynth_t *ref,  char *read_id, uint64_t len_raw_signal){
    // Output of results
    //uint64_t start_event_idx =  db->qstart[i];
    //uint64_t end_event_idx =  db->qend[i];

    ASSERT(start_event_idx>=0 && start_event_idx<=et.n);
    ASSERT(end_event_idx>=0 && end_event_idx<=et.n);

    //fprintf(stderr,"start_event_idx: %ld, end_event_idx: %ld\n", start_event_idx, end_event_idx);


    uint64_t start_raw_idx = et.event[start_event_idx].start; //inclusive
    //fprintf(stderr,"start_raw_idx: %ld\n", start_raw_idx);

    uint64_t end_raw_idx = et.event[end_event_idx].start + et.event[end_event_idx].length; //exclusive
    //fprintf(stderr,"end_raw_idx: %ld\n", et.event[end_event_idx].start);

    uint64_t query_size =  end_event_idx-start_event_idx;
    float block_len = aln.pos_end - aln.pos_st;
    float residue = block_len - aln.score*block_len/(query_size) ;

    // if(aln.score>70){
    //     continue;
    // }

    kstring_t str;
    kstring_t *sp = &str;
    str_init(sp, sizeof(char)*10000);


    sprintf_append(sp,"%s\t",read_id); // read id name
    sprintf_append(sp,"%ld\t%ld\t%ld\t", len_raw_signal, start_raw_idx, end_raw_idx); // Raw signal length, start and end
    sprintf_append(sp,"%c\t",aln.d); // Direction
    sprintf_append(sp,"%s\t",ref->ref_names[aln.rid]); // reference sequence name
    sprintf_append(sp,"%d\t",ref->ref_seq_lengths[aln.rid]); // reference sequence length


    sprintf_append(sp,"%d\t",aln.pos_st); // Reference start
    sprintf_append(sp,"%d\t",aln.pos_end); // Reference end
    sprintf_append(sp,"%d\t",(int)round(residue)); // Number of residues //todo check this
    sprintf_append(sp,"%d\t",(int)round(block_len)); //  Alignment block length //todo check this
    sprintf_append(sp,"%d\t",aln.mapq); // Mapq
    sprintf_append(sp,"tp:A:P\t");
    sprintf_append(sp,"d1:f:%.2f\t",aln.score); // distance of the best match
    sprintf_append(sp,"d2:f:%.2f\t",aln.score2); // distance of the second best matcj
    sprintf_append(sp,"en:i:%d\n",query_size); // number of events in the mapped segment

    return sp->s;
}

aln_t map(refsynth_t *ref, float *raw, int64_t nsample, int polyend, char *read_id, char **sp, sigfish_opt_t opt, int *nevents){
    ASSERT(ref != NULL);
    ASSERT(raw != NULL);
    ASSERT(nsample > 0);
    ASSERT(nsample-polyend >= opt.query_size_sig);
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
        ASSERT((uint64_t)start_idx < et.n);
        end_idx = start_idx + opt.query_size_events;

        if (start_idx + 25 > et.n ){
            fprintf(stderr,"WARNING: not enough events to map - a weird read (<25 events in %ld samples)\n",nsample-polyend);
            start_idx = 0;end_idx = 0;
        } else if(end_idx > et.n){
            fprintf(stderr,"WARNING: Only %ld events in %ld samples\n",et.n-start_idx,nsample-polyend);
            end_idx = et.n-1; //TODO: this is a hack, investigate why this happens as et.n is supposed to be inclusive
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
                update_aln(aln, min_score, j, min_pos-(qlen-1)*rlen, '+', cost, qlen, rlen, 0);
            }
            free(cost);
        }

        free(query);
        update_best_aln(&best_aln, aln, ref);
        free(aln);

        if(best_aln.pos_st >= 0 && sp!=NULL){
            *sp=sprintf_aln(start_idx, end_idx, et, best_aln,  ref, read_id, nsample);
        } else {
            *sp=NULL;
        }
        *nevents = qlen;
    }
    free(et.event);

    return best_aln;
}

// #define SIGFISH_DTW_CUTOFF 70


void test1(sigfish_rstate_t *r, sigfish_state_t *state, int channel, enum sigfish_status *status, int i){
    if (r->len_raw_signal < 30){
        state->status[channel] = status[i] = SIGFISH_MORE;
    } else {
        float sum = 0;
        for(int j=0;j<30;j++) sum += r->raw_signal[j];
        sum /= 30;

        if((int)sum % 2 == 0){
            state->status[channel] = status[i] = SIGFISH_REJECT;
        } else {
            state->status[channel] = status[i] = SIGFISH_CONT;
        }
        fprintf(stderr,"channel: %d, read %d, sum: %f, status: %d\n",channel,r->read_number,sum,status[i]);
    }
}

void test2(sigfish_rstate_t *r, sigfish_state_t *state, int channel, enum sigfish_status *status, int i){

    jnnv3_aparam_t param = JNNV3_R9_ADAPTOR; // may impact peformance
    if (state->opt.pore == OPT_PORE_RNA004) {
        jnnv3_aparam_t atmp = JNNV3_RNA004_ADAPTOR;
        param = atmp;
    }

    int chunk_size = param.chunk_size;
    uint64_t sigfish_min_samples = chunk_size*(param.start_chunks+1);

    //if too short to start detecting adaptor
    if (r->len_raw_signal < sigfish_min_samples){
        state->status[channel] = status[i] = SIGFISH_MORE;
    } else {

        //detect adaptor
        float sum = 0;
        for(int j = chunk_size*param.start_chunks; j < (int)sigfish_min_samples; j++){
            sum += r->raw_signal[j];
        }
        sum /= chunk_size;

        if((int)sum % 2 == 0){
            state->status[channel] = status[i] = SIGFISH_REJECT;
        } else {
            state->status[channel] = status[i] = SIGFISH_CONT;
        }

        //detect polya

        //dtw

    }
}

 int debug = 0;

void decide(sigfish_rstate_t *r, sigfish_state_t *state, int channel, enum sigfish_status *status, int i){

    if (debug == 0) {
        jnnv3_aparam_t param = JNNV3_R9_ADAPTOR; // may impact peformance
        if (state->opt.pore == OPT_PORE_RNA004) {
            jnnv3_aparam_t atmp = JNNV3_RNA004_ADAPTOR;
            param = atmp;
        }
        uint64_t sigfish_min_samples = param.chunk_size*(param.start_chunks+1);

        state->status[channel] = status[i] = SIGFISH_MORE;
        //if too short to start detecting adaptor
        if (r->len_raw_signal >= sigfish_min_samples){

            float *sig_store = r->raw_signal;
            int sig_store_i = r->len_raw_signal;

            int cur_chunk_st = r->cur_chunk_st;
            float *chunk = &sig_store[cur_chunk_st];
            int current_chunk_size = sig_store_i-cur_chunk_st;

            jnnv3_astate_t *s = state->s[channel];
            jnnv3_pstate_t *t = state->t[channel];

            jnnv3_pparam_t pparam = JNNV3_R9_POLYA;
            if (state->opt.pore == OPT_PORE_RNA004) {
                jnnv3_pparam_t ptmp = JNNV3_RNA004_POLYA;
                pparam = ptmp;
            }

            if (s->top == 0){ //enough chunks arrived
                LOG_TRACE("%s","Enough chunks, start to detect adaptor");
                jnnv3_acalc_param(s, param, sig_store, sig_store_i);
                LOG_TRACE("top %f",s->top);
                chunk = sig_store;
                current_chunk_size = sig_store_i;
            }

            if (!s->adapter_found){
                jnnv3_acore(s, param, chunk, current_chunk_size);
                if (s->adapter_found){
                    jnn_pair_t p = s->segs[0];
                    LOG_TRACE("Adapter found at %ld,%ld. sigstore size %d",p.x,p.y,sig_store_i);
                    jnnv3_pcalc_param(t, p, pparam, sig_store, sig_store_i);
                    chunk = &sig_store[p.y];
                    current_chunk_size = sig_store_i-p.y;

                } else {
                    LOG_TRACE("%s","Adapter not found, continue to detect adaptor");
                }
            }

            if(s->adapter_found && !t->polya_found){
                jnnv3_pcore(t, pparam,chunk,current_chunk_size);
            }

            if(t->polya_found){
                ASSERT(s->adapter_found == 1);
                ASSERT(t->seg_i > 0);
                ASSERT(s->seg_i > 0);
                jnn_pair_t polya = t->segs[0];
                jnn_pair_t adapt = s->segs[0];
                int st = polya.y+adapt.y-1;
                int leftover = sig_store_i - st;
                if(leftover >= state->opt.query_size_sig){
                    //fprintf(stderr,"leftover: %d, running DTW\n", leftover);
                    char *read_id;
                    char tmp[100];
                    if(r->read_id){
                        read_id = r->read_id;
                    } else {
                        sprintf(tmp, "read_%d_channel_%d", r->read_number, channel+1);
                        read_id = tmp;
                    }
                    int nevents = state->opt.query_size_events;
                    char **sp = state->debug ? &(state->debug[i]) : NULL;
                    aln_t best_aln=map(state->ref, sig_store, sig_store_i, st, read_id, sp, state->opt, &nevents);
                    ASSERT(nevents>0)
                    //if(state->debug[i])fprintf(stderr,"%s",state->debug[i]);
                    if(best_aln.score < state->opt.dtw_cutoff * (nevents) / state->opt.query_size_events){
                        state->status[channel] = status[i] = SIGFISH_REJECT;
                    } else {
                        state->status[channel] = status[i] = SIGFISH_CONT;
                    }
                } else {
                    //fprintf(stderr,"leftover: %d, waiting for more\n", leftover);
                }
            }
        }

    }
    else if(debug==1){
        test1(r,state,channel,status,i);
    } else if (debug==2){
        test2(r,state,channel,status,i);
    }

}


void process_sigfish_single(sigfish_state_t *state, sigfish_read_t *read_batch, int i){
    int channel = read_batch[i].channel-1;
    ASSERT(channel>=0 && channel < state->num_channels);

    //populate
    sigfish_rstate_t *r = &state->reads[channel];
    if (r->read_number == read_batch[i].read_number){ //same read number
        if(r->c_raw_signal < r->len_raw_signal + read_batch[i].len_raw_signal){
            r->c_raw_signal = r->len_raw_signal + read_batch[i].len_raw_signal;
            r->raw_signal = (float *)realloc(r->raw_signal, r->c_raw_signal*sizeof(float));
            MALLOC_CHK(r->raw_signal);
        }
        memcpy(r->raw_signal+r->len_raw_signal, read_batch[i].raw_signal, read_batch[i].len_raw_signal*sizeof(float));
        r->cur_chunk_st = r->len_raw_signal;
        r->len_raw_signal += read_batch[i].len_raw_signal;
        LOG_TRACE("same read %d len %ld",r->read_number,r->len_raw_signal);
        if(r->read_id){
            ASSERT(strcmp(r->read_id, read_batch[i].read_id)==0);
        }
    } else { //new read number
        r->len_raw_signal = read_batch[i].len_raw_signal;
        r->read_number = read_batch[i].read_number;
        state->status[channel] = state->status_ret[i] = 0;
        LOG_TRACE("new read %d len %ld",r->read_number,r->len_raw_signal);
        if(read_batch[i].read_id){
            r->read_id=realloc(r->read_id, strlen(read_batch[i].read_id)+1);
            MALLOC_CHK(r->read_id);
            strcpy(r->read_id,read_batch[i].read_id);
        }
        jnnv3_astate_t *s = state->s[channel];
        jnnv3_pstate_t *t = state->t[channel];

        jnnv3_aparam_t param = JNNV3_R9_ADAPTOR;
        jnnv3_pparam_t pparam = JNNV3_R9_POLYA;
        if (state->opt.pore == OPT_PORE_RNA004) {
            jnnv3_aparam_t atmp = JNNV3_RNA004_ADAPTOR;
            param = atmp;
            jnnv3_pparam_t ptmp = JNNV3_RNA004_POLYA;
            pparam = ptmp;
        }


        reset_jnnv3_astate(s,param);
        reset_jnnv3_pstate(t,pparam);
        if(r->c_raw_signal < read_batch[i].len_raw_signal){
            r->c_raw_signal = read_batch[i].len_raw_signal;
            r->raw_signal = (float *)realloc(r->raw_signal, r->c_raw_signal*sizeof(float));
            MALLOC_CHK(r->raw_signal);
        }
        memcpy(r->raw_signal, read_batch[i].raw_signal, read_batch[i].len_raw_signal*sizeof(float));
        r->cur_chunk_st = 0;
    }


    //process
    decide(r, state, channel, state->status_ret, i);
}


enum sigfish_status *process_sigfish(sigfish_state_t *state, sigfish_read_t *read_batch, int batch_size){

    double start = realtime();

    enum sigfish_status *status = (enum sigfish_status *)calloc(state->num_channels,sizeof(enum sigfish_status));
    MALLOC_CHK(status);
    state->status_ret = status;
    state->n_rec = batch_size;
    // for(int i=0;i<batch_size;i++){
    //     process_sigfish_single(state, read_batch, i);
    // }
    work_rt(state, read_batch,process_sigfish_single);

    if(state->debug){

        for(int i=0;i<batch_size;i++){
            if(state->status_ret[i]!=SIGFISH_MORE && state->debug[i]){
                fprintf(state->debug_paf,"%s",state->debug[i]);
                free(state->debug[i]);
            }
        }
    }

    double end = realtime();
    VERBOSE("%d queries processed in %f seconds",batch_size,end-start);

    return status;
}