
/* @file genref.c
**

** @@
******************************************************************************/


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>

#include "sigfish.h"
#include "error.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


model_t *gen_ref(char *genome, model_t *pore_model){

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(genome, "r");
    F_CHK(fp,genome);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    while ((l = kseq_read(seq)) >= 0) {


    }

    kseq_destroy(seq);
    gzclose(fp);

}
