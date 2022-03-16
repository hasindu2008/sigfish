/* @file sigfish.h
**
******************************************************************************/

#ifndef SIGFISH_H
#define SIGFISH_H

#include <stdint.h>
#include "slow5/slow5.h"

#define SIGFISH_VERSION "0.1.0"

//model types
#define MODEL_TYPE_NUCLEOTIDE 1
#define MODEL_TYPE_METH 2

#define MAX_KMER_SIZE 6 //maximum k-mer size
#define MAX_NUM_KMER 4096   //maximum number of k-mers in nucleotide model
#define MAX_NUM_KMER_METH 15625 //maximum number of k-mers in methylated model

//default model IDs
#define MODEL_ID_DNA_NUCLEOTIDE 1
#define MODEL_ID_RNA_NUCLEOTIDE 2

/*******************************************************
 * flags related to the user specified options (opt_t) *
 *******************************************************/

#define SIGFISH_RNA 0x001 //if RNA or not
#define SIGFISH_DTW 0x002 //if dtw-std
#define SIGFISH_INV 0x004 //if set, reverse reference events instead of query events
#define SIGFISH_SEC 0x008 //if secondaries are printed
#define SIGFISH_REF 0x010 //map to the whole reference
#define SIGFISH_END 0x020 //map the end of the query

#define SECONDARY_CAP 5 //maximum number of secondary events to print

#define WORK_STEAL 1 //simple work stealing enabled or not (no work stealing mean no load balancing)
#define STEAL_THRESH 1 //stealing threshold

/* a single signal-space event : adapted from taken from scrappie */
typedef struct {
    uint64_t start;
    float length; //todo : cant be made int?
    float mean;
    float stdv;
    //int32_t pos;   //todo : always -1 can be removed
    //int32_t state; //todo : always -1 can be removed
} event_t;

/* event table : adapted from scrappie */
typedef struct {
    size_t n;     //todo : int32_t not enough?
    size_t start; //todo : always 0?
    size_t end;   //todo : always equal to n?
    event_t* event;
} event_table;

/* k-mer model */
typedef struct {
    float level_mean;
    float level_stdv;

#ifdef CACHED_LOG
    float level_log_stdv;     //pre-calculated for efficiency
#endif

#ifdef LOAD_SD_MEANSSTDV
    //float sd_mean;
    //float sd_stdv;
    //float weight;
#endif
} model_t;

typedef struct {
    int32_t num_ref;
    char **ref_names;
    int32_t *ref_lengths;
    int32_t *ref_seq_lengths;

    float **forward;
    float **reverse;
} refsynth_t;


/* scaling parameters for the signal : taken from nanopolish */
typedef struct {
    // direct parameters that must be set
    float scale;
    float shift;
    //float drift; = 0 always?
    float var; // set later when calibrating
    //float scale_sd;
    //float var_sd;

#ifdef CACHED_LOG
    float log_var;    // derived parameters that are cached for efficiency
#endif
    //float scaled_var;
    //float log_scaled_var;
} scalings_t;


/* user specified options */
typedef struct {
    const char* model_file;     //name of the k-mer model file
    const char* meth_model_file;//name of the methylation model file
    uint32_t flag;              //flags
    int32_t batch_size;         //max reads loaded at once: K
    int64_t batch_size_bytes;   //max bytes loaded at once: B

    int32_t num_thread; //t
    int8_t verbosity;
    int32_t debug_break;

    char *region_str; //the region string in format chr:start-end

    int32_t prefix_size;
    int32_t query_size;

} opt_t;

typedef struct {
    int32_t rid;
    int32_t pos;
    float score;
    float score2;
    char d;
    uint8_t mapq;
} aln_t;

/* a batch of read data (dynamic data based on the reads) */
typedef struct {

    int32_t n_rec;
    int32_t capacity_rec;

    char **mem_records;
    size_t *mem_bytes;

    slow5_rec_t **slow5_rec;

    float **current_signal;

    //event table
    event_table* et;

    //scaling
    scalings_t* scalings;

    //results
    aln_t* aln;

    //stats
    int64_t sum_bytes;
    int64_t total_reads; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)



} db_t;



/* core data structure (mostly static data throughout the program lifetime) */
typedef struct {

    //multi region related
    char **reg_list; //the list of regions
    int64_t reg_n;   //number of regions in list
    int64_t reg_i;   //current region being processed

    //clipping coordinates
    int32_t clip_start;
    int32_t clip_end;

    //slow5
    slow5_file_t *sf;

    // models
    model_t* model; //dna model
    model_t* cpgmodel; //cpg model
    uint32_t kmer_size;

    // options
    opt_t opt;

    //realtime0
    double realtime0;

    double load_db_time;
    double process_db_time;
    double output_time;

    //stats //set by output_db
    int64_t sum_bytes;
    int64_t total_reads; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)

    refsynth_t *ref;

} core_t;


/* argument wrapper for the multithreaded framework used for data processing */
typedef struct {
    core_t* core;
    db_t* db;
    int32_t starti;
    int32_t endi;
    void (*func)(core_t*,db_t*,int);
    int32_t thread_index;
#ifdef WORK_STEAL
    void *all_pthread_args;
#endif
#ifdef HAVE_CUDA
    int32_t *ultra_long_reads; //reads that are assigned to the CPU due to the unsuitability to process on the GPU
    double ret1;    //return value
#endif
} pthread_arg_t;

/* return status by the load_db - used for termination when all the data is processed */
typedef struct {
    int32_t num_reads;
    int64_t num_bytes;
} ret_status_t;

/******************************************
 * function prototype for major functions *
 ******************************************/

/* initialise user specified options */
void init_opt(opt_t* opt);

/* initialise the core data structure */
core_t* init_core(const char* fastafile, char *slow5file, opt_t opt, double realtime0);

/* initialise a data batch */
db_t* init_db(core_t* core);

/* load a data batch from disk */
ret_status_t load_db(core_t* dg, db_t* db);

void work_per_single_read(core_t* core,db_t* db, int32_t i);
/* process all reads in the given batch db */
void work_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));

/* process a data batch */
void process_db(core_t* core, db_t* db);

/* align a single read specified by index i*/
void process_single(core_t* core, db_t* db, int32_t i);

/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db);

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(db_t* db);

/* completely free a data batch */
void free_db(db_t* db);

/* free the core data structure */
void free_core(core_t* core,opt_t opt);

#endif
