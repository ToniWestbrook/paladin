#ifndef ALIGN_H_
#define ALIGN_H_

#include "bwa.h"
#include "bwamem.h"
#include "kseq.h"

KSEQ_DECLARE(gzFile)

typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;


static void * process(void *shared, int step, void *_data);
static void update_a(mem_opt_t *opt, const mem_opt_t *opt0);

// 'align' command entry point
int command_align(int argc, char *argv[]);

int renderAlignUsage(const mem_opt_t * passOptions);

// CLEAN
void *kopen(const char *fn, int *_fd);
int kclose(void *a);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);


#endif /* ALIGN_H_ */
