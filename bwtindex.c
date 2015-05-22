/* Contact: Toni Westbrook <anthonyw@wildcats.unh.edu> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwt.h"
#include "utils.h"

#ifdef _DIVBWT
    #include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
    #include "malloc_wrap.h"
#endif

int is_bwt(ubyte_t *T, int n);

// Pack the given byte value into the BWT (pre-interleaved) at the specified index (4 per 32-bit word)
void packValue(bwt_t * passBWT, int passSeqIdx, bwtint_t passValue) {
	int packShift;

	packShift = (8 * (sizeof(uint32_t) - 1)) - (8 * (passSeqIdx % sizeof(uint32_t)));
	passBWT->bwt[passSeqIdx / sizeof(uint32_t)] |= (passValue << packShift);
}

// Unpack byte from from the BWT (pre-interleaved) at the specified index
ubyte_t unpackValue(bwt_t * passBWT, int passSeqIdx) {
	int packShift;

	packShift = (8 * (sizeof(uint32_t) - 1)) - (8 * (passSeqIdx % sizeof(uint32_t)));
	return (passBWT->bwt[passSeqIdx / sizeof(uint32_t)] >> packShift) & 0xFF;
}

// Obtain sequence length from the packed file
int64_t bwa_seq_len(const char *fn_pac) {
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;

	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return pac_len - 1;
}

// Populate the BWT from a packed file
bwt_t * bwt_pac2bwt(const char *fn_pac) {
	bwt_t *bwt;
	ubyte_t * packedBuf, * unpackedBuf;
	int i, packedSize;
	FILE *fp;

	// Initialize BWT structure
	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	bwt->seq_len = bwa_seq_len(fn_pac);
	bwt->bwt_size = (bwt->seq_len + sizeof(uint32_t) - 1) / sizeof(uint32_t);

	fp = xopen(fn_pac, "rb");

	// Prepare buffers
	packedSize = bwt->seq_len;
	packedBuf = (ubyte_t*)calloc(packedSize, 1);
	unpackedBuf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
	memset(bwt->L2, 0, (VALUE_DOMAIN + 1) * sizeof(bwtint_t));

	err_fread_noeof(packedBuf, 1, packedSize, fp);
	err_fclose(fp);

	// Unpack sequence, record occurrence
	for (i = 0; i < bwt->seq_len; ++i) {
		unpackedBuf[i] = packedBuf[i];
		++bwt->L2[unpackedBuf[i] + 1];
	}
	free(packedBuf);

	// Accumulate lower occurrences
	for (i = 2; i <= VALUE_DOMAIN ; ++i) {
		bwt->L2[i] += bwt->L2[i-1];
	}

	// Burrows-Wheeler Transform
	bwt->primary = is_bwt(unpackedBuf, bwt->seq_len);

	// Pack result into BWT
	bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, sizeof(uint32_t));
	for (i = 0; i < bwt->seq_len; ++i) {
		packValue(bwt, i, unpackedBuf[i]);
	}

	free(unpackedBuf);

	return bwt;
}

// the "pac2bwt" command; IMPORTANT: bwt generated at this step CANNOT be used with BWA. bwtupdate is required! 
int bwa_pac2bwt(int argc, char *argv[]) { 
	bwt_t *bwt;
	int c, use_is = 1;
	while ((c = getopt(argc, argv, "d")) >= 0) {
		switch (c) {
		case 'd': use_is = 0; break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa pac2bwt [-d] <in.pac> <out.bwt>\n");
		return 1;
	}
	bwt = bwt_pac2bwt(argv[optind]);
	bwt_dump_bwt(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

// Interleave occurrence counts into the BWT at specified interval for efficient search
void bwt_bwtupdate_core(bwt_t *bwt) {
	bwtint_t i, k, counts[VALUE_DOMAIN], n_occ;
	uint32_t * bwtBuf;

	// Initialize counts array
	memset(counts, 0, sizeof(counts));

	// Adjust the FM-Index size by the number of interleaved occurrence counts
	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN;
	bwtBuf = (uint32_t*)calloc(bwt->bwt_size, sizeof(uint32_t));

	// Iterate through each packed value in the current FM-Index
	for (i = k = 0; i < bwt->seq_len; ++i) {
		// For each occurrence interval, insert the number of running occurrences
		if (i % OCC_INTERVAL == 0) {
			memcpy(bwtBuf + k, counts, sizeof(bwtint_t) * VALUE_DOMAIN);
			k += sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN;
		}

		// Copy packed value
		if (i % 4 == 0) bwtBuf[k++] = bwt->bwt[i/4];

		// Unpack and record current occurrence
		++counts[unpackValue(bwt, i)];
	}

	// Record last element
	memcpy(bwtBuf + k, counts, sizeof(bwtint_t) * VALUE_DOMAIN);
	xassert(k + sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN == bwt->bwt_size, "inconsistent bwt_size");

	// Update FM-Index
	free(bwt->bwt);
	bwt->bwt = bwtBuf;
}

// the "bwtupdate" command
int bwa_bwtupdate(int argc, char *argv[]) {
	bwt_t *bwt;
	if (argc < 2) {
		fprintf(stderr, "Usage: bwa bwtupdate <the.bwt>\n");
		return 1;
	}
	bwt = bwt_restore_bwt(argv[1]);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(argv[1], bwt);
	bwt_destroy(bwt);
	return 0;
}

// the "bwt2sa" command
int bwa_bwt2sa(int argc, char *argv[]) {
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
		case 'i': sa_intv = atoi(optarg); break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa bwt2sa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_sa(bwt, sa_intv);
	bwt_dump_sa(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

// Primary Index command.  Create protein file, pack, construct BWT, interleave, create SA, repack
int bwa_index(int argc, char *argv[]) { 
	bwt_t *bwt;
	char * prefix, * proName, * pacName, * bwtName, * saName;
	gzFile fp;
	clock_t t;
	int64_t l_pac;

	// Setup filenames
	prefix = malloc(strlen(argv[optind]));
	proName = malloc(strlen(argv[optind]) + 5);
	pacName = malloc(strlen(argv[optind]) + 5);
	bwtName = malloc(strlen(argv[optind]) + 5);
	saName = malloc(strlen(argv[optind]) + 5);

	sprintf(prefix, "%s", argv[optind]);
	sprintf(proName, "%s.pro", argv[optind]);
	sprintf(pacName, "%s.pac", argv[optind]);
	sprintf(bwtName, "%s.bwt", argv[optind]);
	sprintf(saName, "%s.sa", argv[optind]);

	// Create Protein Sequence
	t = clock();
	fprintf(stderr, "[panda_protein] Translating protein sequence...\n");
	//writeIndexProtein(prefix, proName, argv[optind+1]);
	writeIndexMultiProtein(prefix, proName, argv[optind+1]);
	fprintf(stderr, "[bwa_protein] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// Pack FASTA
	fp = xzopen(proName, "r");
	t = clock();
	fprintf(stderr, "[panda_index] Packing protein sequence... ");
	l_pac = bns_fasta2bntseq(fp, prefix, 0);
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	err_gzclose(fp);

	// Construct BWT
	t = clock();
	fprintf(stderr, "[panda_index] Constructing BWT for the packed sequence...\n");
	bwt = bwt_pac2bwt(pacName);
	bwt_dump_bwt(bwtName, bwt);
	bwt_destroy(bwt);
	fprintf(stderr, "[panda_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// Update BWT
	t = clock();
	fprintf(stderr, "[panda_index] Update BWT... ");
	bwt = bwt_restore_bwt(bwtName);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(bwtName, bwt);
	bwt_destroy(bwt);
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// Pack Forward-Only FASTA
	fp = xzopen(proName, "r");
	t = clock();
	fprintf(stderr, "[panda_index] Pack forward-only protein sequence... ");
	l_pac = bns_fasta2bntseq(fp, prefix, 1);
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	err_gzclose(fp);

	// Construct Suffix Array from FM-Index and Occurrences
	t = clock();
	fprintf(stderr, "[panda_index] Construct SA from BWT and Occ... ");
	bwt = bwt_restore_bwt(bwtName);
	bwt_cal_sa(bwt, 32);
	bwt_dump_sa(saName, bwt);
	bwt_destroy(bwt);
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	free(prefix);
	free(proName);
	free(pacName);
	free(bwtName);
	free(saName);

	return 0;
}