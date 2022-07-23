/* @file  dtw_main.c
**
** @@
******************************************************************************/
#define _XOPEN_SOURCE 700
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
    {"slow5", required_argument, 0, 's'},          //0 SLOW5 file //removed can be reused
    {"genome", required_argument, 0, 'g'},         //1 reference genome //removed can be reused
    {"threads", required_argument, 0, 't'},        //2 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //3 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},      //4 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},        //5 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //6
    {"version", no_argument, 0, 'V'},              //7
    {"kmer-model", required_argument, 0, 0},       //8 custom nucleotide k-mer model file
    {"meth-model", required_argument, 0, 0},       //9 custom methylation k-mer model file
    {"output",required_argument, 0, 'o'},          //10 output to a file [stdout]
    {"window",required_argument, 0, 'w'},          //11 the genomic window (region)
    {"rna",no_argument,0,0},                       //12 if RNA (eventalign only)
    {"prefix",required_argument,0,'b'},            //13
    {"query-size",required_argument,0,'q'},        //14
    {"debug-break",required_argument, 0, 0},       //15 break after processing the first batch (used for debugging)
    {"dtw-std", no_argument, 0, 0},                //16 Use standard DTW instead of subsequence dtw
    {"invert", no_argument, 0, 0},                 //17 Reverse the reference instead of query
    {"secondary", required_argument, 0, 0},        //18 Print secondary mappings or not
    {"full-ref", no_argument, 0, 0},               //19 Map to full reference instead of a segment
    {"from-end", no_argument, 0, 0},               //20 Map the end portion of the query
    {0, 0, 0, 0}};


static inline int64_t mm_parse_num(const char* str) //taken from minimap2
{
    double x;
    char* p;
    x = strtod(str, &p);
    if (*p == 'G' || *p == 'g')
        x *= 1e9;
    else if (*p == 'M' || *p == 'm')
        x *= 1e6;
    else if (*p == 'K' || *p == 'k')
        x *= 1e3;
    return (int64_t)(x + .499);
}

static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: sigfish dtw [OPTIONS] genome.fa reads.blow5\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -w STR                     limit processing to a genomic region specified as chr:start-end or a list of regions in a .bed file\n");
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bytes loaded at once [%.1fM]\n",opt.batch_size_bytes/(float)(1000*1000));
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   -o FILE                    output to file [stdout]\n");
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");

    fprintf(fp_help,"\nadvanced options:\n");
    fprintf(fp_help,"   --kmer-model FILE          custom nucleotide k-mer model file (format similar to test/r9-models/r9.4_450bps.nucleotide.6mer.template.model)\n");
    fprintf(fp_help,"   --rna                      the dataset is direct RNA\n");
    fprintf(fp_help,"   -q INT                     the number of events in query signal to align [%d]\n",opt.query_size);
    fprintf(fp_help,"   -p INT                     the number of events to trim at query signal start [%d]\n",opt.prefix_size);
    fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
    fprintf(fp_help,"   --dtw-std                  use DTW standard instead of DTW subsequence\n");
    fprintf(fp_help,"   --invert                   reverse the reference events instead of query\n");
    fprintf(fp_help,"   --secondary STR            print secondary mappings. yes or no [%s]\n",(opt.flag&SIGFISH_SEC)?"yes":"no");
    fprintf(fp_help,"   --full-ref                 map to the full reference\n");
    fprintf(fp_help,"   --from-end                 Map the end portion of the query instead of the beginning\n");
}

//parse yes or no arguments : taken from minimap2
static inline void yes_or_no(opt_t* opt, uint64_t flag, int long_idx,
                             const char* arg,
                             int yes_to_set)
{
    if (yes_to_set) {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag |= flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag &= ~flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    } else {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag &= ~flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag |= flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    }
}


int dtw_main(int argc, char* argv[]) {

    double realtime0 = realtime();

    //signal(SIGSEGV, sig_handler);

    const char* optstring = "p:q:t:B:K:v:o:w:hV";

    int longindex = 0;
    int32_t c = -1;

    char* fastafile = NULL;
    char *slow5file = NULL;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c == 'w') {
            opt.region_str = optarg;
        } else if (c == 'B') {
            opt.batch_size_bytes = mm_parse_num(optarg);
            if(opt.batch_size_bytes<=0){
                ERROR("%s","Maximum number of bytes should be larger than 0.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d", opt.num_thread);
                exit(EXIT_FAILURE);
            }
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum sigfish_log_level_opt)v);
        } else if (c=='V'){
            fprintf(stdout,"sigfish %s\n",SIGFISH_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        } else if (c == 'p'){ // prefix size
            opt.prefix_size = atoi(optarg);
            if (opt.prefix_size < 0) {
                INFO("%s","Autodetect query start.");
            }
        } else if (c == 'q'){ //query size
            opt.query_size = atoi(optarg);
            if (opt.query_size < 0) {
                ERROR("Query size should larger than 0. You entered %d",opt.query_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 0 && longindex == 8) { //custom nucleotide model file
            opt.model_file = optarg;
        } else if (c == 0 && longindex == 9){ //custom methylation k-mer model
            opt.meth_model_file = optarg;
        } else if(c=='o'){
			if (strcmp(optarg, "-") != 0) {
				if (freopen(optarg, "wb", stdout) == NULL) {
					ERROR("failed to write the output to file %s : %s",optarg, strerror(errno));
					exit(EXIT_FAILURE);
				}
			}
        } else if (c == 0 && longindex == 12){ //if RNA
            yes_or_no(&opt, SIGFISH_RNA, longindex, "yes", 1);
        } else if(c == 0 && longindex == 15){ //debug break
            opt.debug_break = atoi(optarg);
        } else if(c == 0 && longindex == 16){ //dtw variant
            yes_or_no(&opt, SIGFISH_DTW, longindex, "yes", 1);
        } else if(c == 0 && longindex == 17){ //reverse the reference instead of query
            yes_or_no(&opt, SIGFISH_INV, longindex, "yes", 1);
        } else if(c == 0 && longindex == 18){ //secondary mappings
            yes_or_no(&opt, SIGFISH_SEC, longindex, optarg, 1);
        } else if(c == 0 && longindex == 19){ //use full reference
            yes_or_no(&opt, SIGFISH_REF, longindex, "yes", 1);
        } else if(c == 0 && longindex == 20){ //map query end
            yes_or_no(&opt, SIGFISH_END, longindex, "yes", 1);
        }
    }

    // No arguments given
    if (argc - optind != 2 || fp_help == stdout) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    fastafile = argv[optind];
    slow5file = argv[optind+1];

    if (!(opt.flag & SIGFISH_RNA)){ //dna
        if(opt.flag & SIGFISH_DTW){
            ERROR("%s","DTW is only available for RNA.");
            exit(EXIT_FAILURE);
        }
        if(opt.flag & SIGFISH_INV){
            ERROR("%s","Inversion is only available for RNA.");
            exit(EXIT_FAILURE);
        }
        if(opt.flag & SIGFISH_REF){
            ERROR("%s","--full-ref is only available for RNA.");
            exit(EXIT_FAILURE);
        }
    }

    if (opt.prefix_size < 0){
        if(!(opt.flag & SIGFISH_RNA)){
            ERROR("%s","DNA does not support auto query start detection.");
            exit(EXIT_FAILURE);
        } else {
            if(opt.flag & SIGFISH_INV){
                ERROR("%s","Inversion is not compatible with auto query start detection.");
                exit(EXIT_FAILURE);
            }
            if(opt.flag & SIGFISH_END){
                ERROR("%s","Mapping from query end is not compatible with auto query start detection.");
                exit(EXIT_FAILURE);
            }
        }
    }

    if (slow5file == NULL || fastafile == NULL || fp_help == stdout) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    //initialise the core data structure
    core_t* core = init_core(fastafile, slow5file, opt, realtime0);

    int32_t counter=0;

    //initialise a databatch
    db_t* db = init_db(core);

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bytes};
    while (status.num_reads >= core->opt.batch_size || status.num_bytes>=core->opt.batch_size_bytes) {

        //load a databatch
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bytes/(1000.0*1000.0));

        //process a databatch
        process_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bytes/(1000.0*1000.0));

        //output print
        output_db(core, db);

        //free temporary
        free_db_tmp(db);

        if(opt.debug_break==counter){
            break;
        }
        counter++;
    }

    //free the databatch
    free_db(db);

    fprintf(stderr, "[%s] total entries: %ld\tprefix fail: %ld\tignored: %ld\ttoo short: %ld", __func__,(long)core->total_reads, (long)core->prefix_fail, (long)core->ignored, (long)core->too_short);
    fprintf(stderr,"\n[%s] total bytes: %.1f M",__func__,core->sum_bytes/(float)(1000*1000));

    fprintf(stderr, "\n[%s] Data loading time: %.3f sec", __func__,core->load_db_time);
    fprintf(stderr, "\n[%s] Data processing time: %.3f sec", __func__,core->process_db_time);
    fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);

    fprintf(stderr,"\n");

    //free the core data structure
    free_core(core,opt);


    return 0;
}
