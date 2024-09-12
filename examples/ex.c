#include <sigfish.h>

#define CHANNELS 10
#define CHUNK_SIZE 1200 // todo: this variable may be depended on pore
#define ROUNDS 10

int main(){

    sigfish_opt_t opt;
    opt.num_thread = 1;
    opt.debug_paf = NULL;
    opt.no_full_ref = 0;
    opt.dtw_cutoff = 70.0;
    opt.query_size_events = 250;
    opt.query_size_sig = 6000;
    opt.pore = 0;
    opt.samples_per_event = 24;

    sigfish_state_t *state = init_sigfish(NULL, CHANNELS, opt);
    sigfish_read_t reads[CHANNELS];

    for(int r=0; r<ROUNDS; r++){
        printf("round %d\n", r);
        for (int i=0; i<CHANNELS; i++){
            reads[i].channel = i+1;
            reads[i].read_number = 0;
            reads[i].len_raw_signal = CHUNK_SIZE;
            reads[i].raw_signal = (float*)malloc(sizeof(float)*CHUNK_SIZE);
            reads[i].read_id = NULL;
            for (int j=0; j<CHUNK_SIZE; j++){
                reads[i].raw_signal[j] = j+i+r;
            }

        }
        enum sigfish_status *status = process_sigfish(state, reads, CHANNELS);
        for(int i=0;i<CHANNELS;i++){
            printf("channel %d: %d\n", i+1, status[i]);
        }
        putc('\n', stdout);
        free(status);

        for(int i=0; i<CHANNELS; i++){
            free(reads[i].raw_signal);
        }

    }

    free_sigfish(state);

    return 0;
}