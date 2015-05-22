#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwt.h"
#include "utils.h"
#include "kseq.h"

KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_MAP_INIT_STR(str, int)

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

//#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
//#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
//#define kv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
/*
uint32_t * getOccInterval(bwt_t * passBWT, bwtint_t passIndex) {
	int numIntervals;

	numIntervals = passIndex / 128;

	// Each interval stores 64-bit occurrences for each of the values + 128 packed values between
	return passBWT->bwt + (numIntervals * 2 * VALUE_DOMAIN) + (numIntervals * (128 / 4));
}

void packValue(bwt_t * passBWT, int passSeqIdx, bwtint_t passValue) {
	int packShift;

	packShift = (8 * (sizeof(uint32_t) - 1)) - (8 * (passSeqIdx % sizeof(uint32_t)));
	passBWT->bwt[passSeqIdx / sizeof(uint32_t)] |= (passValue << packShift);
}

ubyte_t unpackValue(bwt_t * passBWT, int passSeqIdx) {
	int packShift;

	packShift = (8 * (sizeof(uint32_t) - 1)) - (8 * (passSeqIdx % sizeof(uint32_t)));
	return (passBWT->bwt[passSeqIdx / sizeof(uint32_t)] >> packShift) & 0xFF;
}

ubyte_t unpackBWTValue(bwt_t * passBWT, int passSeqIdx) {
	int packShift;
	uint32_t * bwtValueStart;

	packShift = (8 * (sizeof(uint32_t) - 1)) - (8 * (passSeqIdx % sizeof(uint32_t)));

	bwtValueStart = getOccInterval(passBWT, passSeqIdx);
	bwtValueStart += sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN;

	return (bwtValueStart[(passSeqIdx % 128) / sizeof(uint32_t)] >> packShift) & 0xFF;
}

// DONE
static inline int TEST__occ_aux(uint64_t y, int c) {
	unsigned char * bytePtr;
	int count, byteIdx;

	bytePtr = (unsigned char *) &y;
	count = 0;

	for (byteIdx = 0 ; byteIdx < sizeof(y) ; byteIdx++) {
		if (bytePtr[byteIdx] == c) count++;
	}

	return count;
}

//DONE
uint8_t * TEST_add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q) {
	bntann1_t *p;
	int i, lasts;
	if (bns->n_seqs == *m_seqs) {
		*m_seqs <<= 1;
		bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
	}
	p = bns->anns + bns->n_seqs;
	p->name = strdup((char*)seq->name.s);
	p->anno = seq->comment.l > 0? strdup((char*)seq->comment.s) : strdup("(null)");
	p->gi = 0; p->len = seq->seq.l;
	p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
	p->n_ambs = 0;

	// Iterate through each AA in the sequence
	for (i = lasts = 0; i < seq->seq.l; ++i) {
		// Treat IUPAC as 'A' indexed ASCII equivalent
		int c = (int)(seq->seq.s[i] - 'A');

		// X considered unknown
		if (c == 23) {
			// Contiguous X
			if (lasts == seq->seq.s[i]) {
				++(*q)->len;
			}
			else {
				if (bns->n_holes == *m_holes) {
					(*m_holes) <<= 1;
					bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
				}

				*q = bns->ambs + bns->n_holes;
				(*q)->len = 1;
				(*q)->offset = p->offset + i;
				(*q)->amb = seq->seq.s[i];
				++p->n_ambs;
				++bns->n_holes;
			}
		}

		// Store value to check for contiguous runs
		lasts = seq->seq.s[i];

		// Expand buffer if necessary
		if (bns->l_pac == *m_pac) { // double the pac size
			*m_pac <<= 1;
			pac = realloc(pac, *m_pac);
			memset(pac + bns->l_pac, 0, (*m_pac - bns->l_pac));
		}

		// Pack and set value
		if (c == 23) c = lrand48()&3; // Fix this
		pac[bns->l_pac] = c;

		++bns->l_pac;
	}

	++bns->n_seqs;
	return pac;
}

// DONE
bwtint_t TEST_bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c) {
	bwtint_t n;
	uint32_t *p, *end;

	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (bwtint_t)(-1)) return 0;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// Retrieve Occurrences at k/OCC_INTERVAL, then advance to first BWT cell
	n = ((bwtint_t*)(p = getOccInterval(bwt, k)))[c];
	p += sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN;

	// calculate Occ up to the last k/32
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += TEST__occ_aux((uint64_t)p[0]<<32 | p[1], c);

	// calculate Occ
	n += TEST__occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31; // corrected for the masked bits

	return n;
}

// DONE
static inline bwtint_t TEST_bwt_invPsi(const bwt_t *bwt, bwtint_t k) {
	bwtint_t x = k - (k > bwt->primary);
	x = unpackBWTValue(bwt, x);
	x = bwt->L2[x] + TEST_bwt_occ(bwt, k, x);

	return k == bwt->primary? 0 : x;
}

// DONE
int64_t TEST_bwa_seq_len(const char *fn_pac) {
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

// DONE
void TEST_bwt_dump_bwt(const char *fn, const bwt_t *bwt) {
	FILE *fp;

	fp = xopen(fn, "wb");

	// Format: <primary><L2[1::]><bwt>
	// Sizes: 8, Domain * 8, BWTSize * 4
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), VALUE_DOMAIN, fp);
	err_fwrite(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp);
	err_fflush(fp);
	err_fclose(fp);
}

// DONE
void TEST_bwt_dump_sa(const char *fn, const bwt_t *bwt) {
	FILE *fp;

	fp = xopen(fn, "wb");

	// Format: <primary><L2[1::]><SA Interval><Seq Length><SA>
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), VALUE_DOMAIN, fp);
	err_fwrite(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fwrite(&bwt->seq_len, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	err_fflush(fp);
	err_fclose(fp);
}

//DONE
bwt_t * TEST_bwt_restore_bwt(const char *fn) {
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = xopen(fn, "rb");
	err_fseek(fp, 0, SEEK_END);

	// Populate BWT Size and allocate BWT
	bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * (VALUE_DOMAIN + 1)) / sizeof(uint32_t);
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, sizeof(uint32_t));

	// Populate Primary, L2[1::], and BWT
	err_fseek(fp, 0, SEEK_SET);
	err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), VALUE_DOMAIN, fp);
	fread_fix(fp, bwt->bwt_size * sizeof(uint32_t), bwt->bwt);

	// Sequence length in end of L2 array
	bwt->seq_len = bwt->L2[VALUE_DOMAIN];

	// Generate CNT table
	bwt_gen_cnt_table(bwt);

	printf("restore info - seqlen (%d) bwtlen (%d)\n", bwt->seq_len, bwt->bwt_size);
	int lcv;
	for (lcv = 0 ; lcv < bwt->bwt_size ; lcv++ ) {
		printf("Restored[%d]=%x\n", lcv, bwt->bwt[lcv]);
	}

	err_fclose(fp);

	return bwt;
}

//DONE
void TEST_bwt_bwtupdate_core(bwt_t *bwt) {
	bwtint_t i, k, counts[VALUE_DOMAIN], n_occ;
	uint32_t * bwtBuf;

	// Initialize counts array
	memset(counts, 0, sizeof(counts));

	int lcv;
	for (lcv = 0 ; lcv < bwt->seq_len ; lcv++) {
		printf("macro reports[%d] = %d\n", lcv, unpackValue(bwt, lcv));
	}

	// Adjust the FM-Index size by the number of interleaved occurrence counts
	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN;
	printf("bwt_size %d\n", bwt->bwt_size);
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

		for (lcv = 0 ; lcv < VALUE_DOMAIN ; lcv++) {
			printf("update[%d]=%d\n", lcv, counts[lcv]);
		}
	}

	// Record last element
	memcpy(bwtBuf + k, counts, sizeof(bwtint_t) * VALUE_DOMAIN);
	xassert(k + sizeof(bwtint_t) == bwt->bwt_size, "inconsistent bwt_size");

	// Update FM-Index
	free(bwt->bwt); bwt->bwt = bwtBuf;

	int test = bwt->bwt_size * 4;
	for (lcv = 0 ; lcv < test ; lcv++) {
		printf("FMIndex[%d]=%d\n", lcv, unpackValue(bwt, lcv));

	}

}

// DONE
void TEST_bwt_cal_sa(bwt_t *bwt, int intv) {
	bwtint_t isa, sa, i; // S(isa) = sa
	int intv_round = intv;

	// Sanity checking
	kv_roundup32(intv_round);
	xassert(intv_round == intv, "SA sample interval is not a power of 2.");
	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	// Initialize Suffix Array
	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));

	// Calculate Suffix Array Value
	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % intv == 0) bwt->sa[isa/intv] = sa;
		--sa;
		isa = TEST_bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;

	// SA[0] is currently set to sequence length - set to maximum value
	bwt->sa[0] = (bwtint_t)-1;
}

//Done
int64_t TEST_bns_fasta2bntseq(gzFile fp_fa, const char *prefix, int for_only) {
	extern void seq_reverse(int len, ubyte_t *seq, int is_comp); // in bwaseqio.c
	kseq_t *seq;
	char name[1024];
	bntseq_t *bns;
	uint8_t *pac = 0;
	int32_t m_seqs, m_holes;
	int64_t ret = -1, m_pac, l;
	bntamb1_t *q;
	FILE *fp;

	// Initialization of sequence related structures
	seq = kseq_init(fp_fa);
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	bns->seed = 11;
	srand48(bns->seed);
	m_seqs = m_holes = 8; m_pac = 0x10000;
	bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
	bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
	q = bns->ambs;

	// Allocate packed buffer to size of sequence (for now - to be converted to 40-bit word)
	pac = calloc(m_pac, 1);

	// Open packed file for writing
	strcpy(name, prefix); strcat(name, ".pac");
	fp = xopen(name, "wb");

	// Read sequences into pac buffer
	while (kseq_read(seq) >= 0) pac = TEST_add1(seq, bns, pac, &m_pac, &m_seqs, &m_holes, &q);

	// Add reverse complement if not forward only
	if (!for_only) {
		m_pac = bns->l_pac * 2;
		pac = realloc(pac, m_pac);
		memset(pac + bns->l_pac, 0, m_pac - bns->l_pac);

		for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac) {
			pac[bns->l_pac] = VALUE_DOMAIN - 1 - pac[l];
		}
	}

	ret = bns->l_pac;

	// Write pac file
	ubyte_t ct;
	err_fwrite(pac, 1, bns->l_pac, fp);

	// Pad to nearest 40-bit word
	ct = 0;
	err_fwrite(&ct, 1, 1, fp);
	/*
	// the following codes make the pac file size always (l_pac/4+1+1)
	if (bns->l_pac % 4 == 0) {
	}
	ct = bns->l_pac % 4;
	*/
/*
	ct = 0;
	err_fwrite(&ct, 1, 1, fp);

	// Close pac file
	err_fflush(fp);
	err_fclose(fp);

	bns_dump(bns, prefix);
	bns_destroy(bns);
	kseq_destroy(seq);
	free(pac);
	return ret;
}
//DONE
bwt_t * TEST_bwt_pac2bwt(const char *fn_pac) {
	bwt_t *bwt;
	ubyte_t * packedBuf, * unpackedBuf;
	int i, packedSize;
	FILE *fp;

	// Initialize BWT structure
	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	bwt->seq_len = TEST_bwa_seq_len(fn_pac);
	bwt->bwt_size = (bwt->seq_len + sizeof(uint32_t) - 1) / sizeof(uint32_t);
	printf("Seq len: %d, BWT SIZE: %d\n", bwt->seq_len, bwt->bwt_size);
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
		printf("Unpacked Value at %d: %d\n", i, unpackedBuf[i]);
		++bwt->L2[unpackedBuf[i] + 1];
	}
	free(packedBuf);

	// Accumulate lower occurrences
	for (i = 2; i <= VALUE_DOMAIN ; ++i) {
		bwt->L2[i] += bwt->L2[i-1];
		printf("Accmulating - Value of %d: %d\n", i, bwt->L2[i]);
	}

	// Burrows-Wheeler Transform
	bwt->primary = is_bwt(unpackedBuf, bwt->seq_len);

	int testIdx;
	for (testIdx = 0 ; testIdx < bwt->seq_len + 1 ; testIdx++) {
		printf("BWT[%d]=%d\n", testIdx, unpackedBuf[testIdx]);
	}

	// Pack result into BWT
	bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, sizeof(uint32_t));
	for (i = 0; i < bwt->seq_len; ++i) {
		packValue(bwt, i, unpackedBuf[i]);
		printf("after:%x\n", bwt->bwt[i>>4]);
	}
	free(unpackedBuf);

	return bwt;
}
//done
int TEST_bwa_index(int argc, char *argv[]) {
	bwt_t *bwt;
	char * prefix, * pacName, * bwtName, * saName;
	gzFile fp;
	clock_t t;
	int64_t l_pac;

	// Setup filenames
	prefix = malloc(strlen(argv[optind]));
	pacName = malloc(strlen(argv[optind]) + 4);
	bwtName = malloc(strlen(argv[optind]) + 4);
	saName = malloc(strlen(argv[optind]) + 4);

	sprintf(prefix, "%s", argv[optind]);
	sprintf(pacName, "%s.pac", argv[optind]);
	sprintf(bwtName, "%s.bwt", argv[optind]);
	sprintf(saName, "%s.sa", argv[optind]);

	// Pack FASTA
	fp = xzopen(prefix, "r");
	t = clock();
	fprintf(stderr, "[bwa_index] Packing amino acid sequence... ");
	l_pac = TEST_bns_fasta2bntseq(fp, prefix, 0);
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	err_gzclose(fp);

	// Construct BWT
	t = clock();
	fprintf(stderr, "[bwa_index] Constructing BWT for the packed sequence...\n");
	bwt = TEST_bwt_pac2bwt(pacName);
	TEST_bwt_dump_bwt(bwtName, bwt);
	bwt_destroy(bwt);
	fprintf(stderr, "[bwa_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// Update BWT
	t = clock();
	fprintf(stderr, "[bwa_index] Update BWT... ");
	bwt = TEST_bwt_restore_bwt(bwtName);
	TEST_bwt_bwtupdate_core(bwt);
	TEST_bwt_dump_bwt(bwtName, bwt);
	bwt_destroy(bwt);
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// Pack Forward-Only FASTA
	fp = xzopen(argv[optind], "r");
	t = clock();
	fprintf(stderr, "[bwa_index] Pack forward-only FASTA... ");
	l_pac = TEST_bns_fasta2bntseq(fp, prefix, 1);
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	err_gzclose(fp);

	// Construct Suffix Array from FM-Index and Occurrences
	t = clock();
	fprintf(stderr, "[bwa_index] Construct SA from BWT and Occ... ");
	bwt = TEST_bwt_restore_bwt(bwtName);
	TEST_bwt_cal_sa(bwt, 32);
	TEST_bwt_dump_sa(saName, bwt);
	bwt_destroy(bwt);
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

}


bwaidx_t *bwa_idx_load_from_disk(const char *hint, int which)
{
	bwaidx_t *idx;
	char *prefix;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	idx = calloc(1, sizeof(bwaidx_t));
	if (which & BWA_IDX_BWT) idx->bwt = bwa_idx_load_bwt(hint);
	if (which & BWA_IDX_BNS) {
		int i, c;
		idx->bns = bns_restore(prefix);
		for (i = c = 0; i < idx->bns->n_seqs; ++i)
			if (idx->bns->anns[i].is_alt) ++c;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d ALT contigs\n", __func__, c);
		if (which & BWA_IDX_PAC) {
			idx->pac = calloc(idx->bns->l_pac+1, 1);
			err_fread_noeof(idx->pac, 1, idx->bns->l_pac+1, idx->bns->fp_pac);
			err_fclose(idx->bns->fp_pac);
			idx->bns->fp_pac = 0;
		}
	}
	free(prefix);
	return idx;
}

*/



