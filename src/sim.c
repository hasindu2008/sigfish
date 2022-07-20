/* @file  eval.c
**
** @@
******************************************************************************/

#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include "sigfish.h"
#include "error.h"
#include "ref.h"
#include "model.h"
#include "misc.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct{
    char **ref_names;
    int32_t *ref_lengths;
    char **ref_seq;
    int num_ref;
} ref_t;


typedef struct{
    double m;
    double s;
    int64_t x;
} nrng_t;

typedef struct {
    double digitisation;
    double sample_rate;
    double bases_per_second;
    double range;
} profile_t;

profile_t minion_prof = {8192, 4000, 450, 1402.882324};
profile_t prom_prof = {2048, 4000, 450, 748.5801};

typedef struct {
    nrng_t *rand_amp;
    nrng_t *rand_time;
    nrng_t *rand_rlen;
    profile_t profile;
} core_sim_t;


static nrng_t* init_nrng(int64_t seed, double mean, double std){
    nrng_t *rng = (nrng_t *)malloc(sizeof(nrng_t));
    rng->m = mean;
    rng->s = std;
    rng->x = seed;
    return rng;
}
static void free_nrng(nrng_t *rng){
    free(rng);
}

static double rng(nrng_t *r){
    int64_t x = r->x;
    int64_t x_new = (16807 * (x % 127773)) - (2836 * (x / 127773));
    if (x_new > 0) x = x_new; else x = x_new + 2147483647;
    r->x = x_new;
    return (double)x/2147483647;
}

static double nrng(nrng_t *r){
    double u = 0.0;
    double t = 0.0;
    while (u == 0.0) u = rng(r);
    while (t == 0.0) t = 2.0 * 3.14159265 * rng(r);
    double x = sqrt(-2.0 * log(u)) * cos(t);
    return ((x * r->s) + r->m);
}

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {"version", no_argument, 0, 'V'},              //2
    {"output",required_argument, 0, 'o'},          //3 output to a file [stdout]
    {0, 0, 0, 0}};

static ref_t *load_ref(const char *genome){

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(genome, "r");
    F_CHK(fp,genome);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);


    ref_t *ref = (ref_t *) malloc(sizeof(ref_t));

    int c = 1;
    ref->ref_lengths = (int32_t *) malloc(sizeof(int32_t));
    ref->ref_names = (char **) malloc(sizeof(char *));
    ref->ref_seq = (char **) malloc(sizeof(char *));

    int i = 0;
    while ((l = kseq_read(seq)) >= 0) {
        assert(l==(int)strlen(seq->seq.s));

        if(i+1 > c){
            c *= 2;
            ref->ref_lengths = (int32_t *) realloc(ref->ref_lengths, c*sizeof(int32_t));
            ref->ref_names = (char **) realloc(ref->ref_names, c*sizeof(char *));
            ref->ref_seq = (char **) realloc(ref->ref_seq, c*sizeof(char *));
        }


        ref->ref_lengths[i] = l;
        ref->ref_names[i] = (char *) malloc(strlen(seq->name.s)+1);
        strcpy(ref->ref_names[i], seq->name.s);
        ref->ref_seq[i] = (char *) malloc((l+1)*sizeof(char));
        strcpy(ref->ref_seq[i], seq->seq.s);

        i++;

    }

    ref->num_ref = i;

    kseq_destroy(seq);
    gzclose(fp);

    return ref;

}

static void free_ref(ref_t *ref){

    for(int i=0;i<ref->num_ref;i++){
        free(ref->ref_names[i]);
        free(ref->ref_seq[i]);
    }

    free(ref->ref_lengths);
    free(ref->ref_names);
    free(ref->ref_seq);
    free(ref);
}



static void set_header_attributes(slow5_file_t *sp){

    slow5_hdr_t *header=sp->header;

    //add a header group attribute called run_id
    if (slow5_hdr_add("run_id", header) < 0){
        fprintf(stderr,"Error adding run_id attribute\n");
        exit(EXIT_FAILURE);
    }
    //add another header group attribute called asic_id
    if (slow5_hdr_add("asic_id", header) < 0){
        fprintf(stderr,"Error adding asic_id attribute\n");
        exit(EXIT_FAILURE);
    }
    //add another header group attribute called asic_id
    if (slow5_hdr_add("exp_start_time", header) < 0){
        fprintf(stderr,"Error adding asic_id attribute\n");
        exit(EXIT_FAILURE);
    }
    //add another header group attribute called flow_cell_id
    if (slow5_hdr_add("flow_cell_id", header) < 0){
        fprintf(stderr,"Error adding flow_cell_id attribute\n");
        exit(EXIT_FAILURE);
    }

    //set the run_id attribute to "run_0" for read group 0
    if (slow5_hdr_set("run_id", "run_0", 0, header) < 0){
        fprintf(stderr,"Error setting run_id attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }
    //set the asic_id attribute to "asic_0" for read group 0
    if (slow5_hdr_set("asic_id", "asic_id_0", 0, header) < 0){
        fprintf(stderr,"Error setting asic_id attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }
    //set the exp_start_time attribute to "2022-07-20T00:00:00Z" for read group 0
    if (slow5_hdr_set("exp_start_time", "2022-07-20T00:00:00Z", 0, header) < 0){
        fprintf(stderr,"Error setting exp_start_time attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }
    //set the flow_cell_id attribute to "FAN00000" for read group 0
    if (slow5_hdr_set("flow_cell_id", "FAN00000", 0, header) < 0){
        fprintf(stderr,"Error setting flow_cell_id attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }

}

static void set_header_aux_fields(slow5_file_t *sp){

    //add auxilliary field: channel number
    if (slow5_aux_add("channel_number", SLOW5_STRING, sp->header) < 0){
        fprintf(stderr,"Error adding channel_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    //add auxilliary field: median_before
    if (slow5_aux_add("median_before", SLOW5_DOUBLE, sp->header) < 0) {
        fprintf(stderr,"Error adding median_before auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    //add auxilliary field: read_number
    if(slow5_aux_add("read_number", SLOW5_INT32_T, sp->header) < 0){
        fprintf(stderr,"Error adding read_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }
    //add auxilliary field: start_mux
    if(slow5_aux_add("start_mux", SLOW5_UINT8_T, sp->header) < 0){
        fprintf(stderr,"Error adding start_mux auxilliary field\n");
        exit(EXIT_FAILURE);
    }
    //add auxilliary field: start_time
    if(slow5_aux_add("start_time", SLOW5_UINT64_T, sp->header) < 0){
        fprintf(stderr,"Error adding start_time auxilliary field\n");
        exit(EXIT_FAILURE);
    }

}

static void set_record_primary_fields(profile_t *profile, slow5_rec_t *slow5_record, char *read_id, double offset, double range, int64_t len_raw_signal, int16_t *raw_signal){

    slow5_record -> read_id = read_id;
    slow5_record-> read_id_len = strlen(slow5_record -> read_id);
    slow5_record -> read_group = 0;
    slow5_record -> digitisation = profile->digitisation;
    slow5_record -> offset = offset;
    slow5_record -> range = range;
    slow5_record -> sampling_rate = profile->sample_rate;
    slow5_record -> len_raw_signal = len_raw_signal;
    slow5_record -> raw_signal = raw_signal;

}

static void set_record_aux_fields(slow5_rec_t *slow5_record, slow5_file_t *sp, double median_before, int32_t read_number, uint64_t start_time){

    const char *channel_number = "0";
    uint8_t start_mux = 0;

    if(slow5_aux_set_string(slow5_record, "channel_number", channel_number, sp->header) < 0){
        fprintf(stderr,"Error setting channel_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "median_before", &median_before, sp->header) < 0){
        fprintf(stderr,"Error setting median_before auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "read_number", &read_number, sp->header) < 0){
        fprintf(stderr,"Error setting read_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "start_mux", &start_mux, sp->header) < 0){
        fprintf(stderr,"Error setting start_mux auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "start_time", &start_time, sp->header) < 0){
        fprintf(stderr,"Error setting start_time auxilliary field\n");
        exit(EXIT_FAILURE);
    }


}


int16_t *gen_sig(profile_t *profile, const char *read, int32_t len, model_t *pore_model, uint32_t kmer_size, double *range, double *offset, int64_t *len_raw_signal){
    int64_t n_kmers = len-kmer_size+1;
    int64_t n=0;
    int sps = round((profile->sample_rate)/(double)(profile->bases_per_second));
    int64_t c = n_kmers * sps;
    int16_t *raw_signal = (int16_t *)malloc(c*sizeof(int16_t));
    *offset = 4;
    *range = profile->range;
    for (int i=0; i< n_kmers; i++){
        uint32_t kmer_rank = get_kmer_rank(read+i, kmer_size);
        float s = pore_model[kmer_rank].level_mean;
        for(int j=0; j<sps; j++){
            if(n==c){
                c *= 2;
                raw_signal = (int16_t *)realloc(raw_signal, c*sizeof(int16_t));
            }
            raw_signal[n] = s*(profile->digitisation)/(*range)-(*offset);
            n++;
        }
    }
    *len_raw_signal = n;
    return raw_signal;
}


static core_sim_t *init_core(){
    core_sim_t *core = (core_sim_t *)malloc(sizeof(core_sim_t));
    core->profile = minion_prof;
    return core;
}

void free_core(core_sim_t *core){
    free(core);
}

int sim_main(int argc, char* argv[]) {

    const char* optstring = "o:hV";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    char *output_file = NULL;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='V'){
            fprintf(stdout,"sigfish %s\n",SIGFISH_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if(c == 'o'){
            output_file=optarg;
        }
    }

    if (argc-optind<1 || output_file==NULL ||  fp_help == stdout) {
        fprintf(fp_help,"Usage: sigfish sim ref.fa\n");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -h                         help\n");
        fprintf(fp_help,"   --version                  print version\n");
        fprintf(fp_help,"   -o FILE                    SLOW5/BLOW5 file to write.\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    char *refname = argv[optind];
    ref_t *ref = load_ref(refname);

    model_t *model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER);
    uint32_t kmer_size = set_model(model, MODEL_ID_DNA_NUCLEOTIDE);

    slow5_file_t *sp = slow5_open(output_file, "w");
    if(sp==NULL){
        fprintf(stderr,"Error opening file!\n");
        exit(EXIT_FAILURE);
    }

    set_header_attributes(sp);
    set_header_aux_fields(sp);

    if(slow5_hdr_write(sp) < 0){
        fprintf(stderr,"Error writing header!\n");
        exit(EXIT_FAILURE);
    }

    int n = 1;
    double median_before = 400;
    int64_t n_samples = 0;
    double range = 0;
    double offset = 0;
    int64_t len_raw_signal =0;

    core_sim_t *core = init_core();

    for(int i=0;i<n;i++){

        slow5_rec_t *slow5_record = slow5_rec_init();
        if(slow5_record == NULL){
            fprintf(stderr,"Could not allocate space for a slow5 record.");
            exit(EXIT_FAILURE);
        }
        int16_t *raw_signal =gen_sig(&core->profile, ref->ref_seq[0], ref->ref_lengths[0], model, kmer_size, &range, &offset, &len_raw_signal);

        char *read_id= (char *)malloc(sizeof(char)*(1000));
        sprintf(read_id,"read_%d",i);
        printf(">%s\n",read_id);
        printf("%s\n",ref->ref_seq[0]);


        set_record_primary_fields(&core->profile, slow5_record, read_id, 6,140, len_raw_signal, raw_signal);
        set_record_aux_fields(slow5_record, sp, median_before, i, n_samples);
        n_samples+=len_raw_signal;


        if (slow5_write(slow5_record, sp) < 0){
            fprintf(stderr,"Error writing record!\n");
            exit(EXIT_FAILURE);
        }

        slow5_rec_free(slow5_record);

        fprintf(stderr,"%d reads done\n",i+1);

    }

    free_core(core);
    free(model);
    slow5_close(sp);
    free_ref(ref);

    return 0;
}
