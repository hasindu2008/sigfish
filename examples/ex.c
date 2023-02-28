#include <sigfish.h>

#define CHANNELS 10
#define CHUNK_SIZE 1200
#define ROUNDS 10

int main(){

    sigfish_state_t *state = init_sigfish(NULL, CHANNELS, 1);
    sigfish_read_t reads[CHANNELS];

    for(int r=0; r<ROUNDS; r++){
        printf("round %d\n", r);
        for (int i=0; i<CHANNELS; i++){
            reads[i].channel = i;
            reads[i].read_number = 0;
            reads[i].len_raw_signal = CHUNK_SIZE;
            reads[i].raw_signal = (float*)malloc(sizeof(float)*CHUNK_SIZE);
            for (int j=0; j<CHUNK_SIZE; j++){
                reads[i].raw_signal[j] = j+i+r;
            }

        }
        enum sigfish_status *status = process_sigfish(state, reads, CHANNELS);
        for(int i=0;i<CHANNELS;i++){
            printf("channel %d: %d\n", i, status[i]);
        }
        putc('\n', stdout);
        free(status);

        for(int i=0; i<10; i++){
            free(reads[i].raw_signal);
        }

    }

    free_sigfish(state);

    return 0;
}