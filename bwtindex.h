#ifndef BWTINDEX_H_
#define BWTINDEX_H_

#include "bwt.h"

#define INDEX_COMPATIBILITY_NONE 0
#define INDEX_COMPATIBLITY_FULL 1
#define INDEX_COMPATBILITY_FUTURE 2

typedef struct {
	int nucleotide;
	int multiFrame;
	int referenceType;
	int version[3];
} IndexHeader;

int command_index(int argc, char *argv[]);
int command_prepare(int argc, char *argv[]);

// Header IO
void writeIndexHeader(FILE * passFilePtr, IndexHeader passHeader);
IndexHeader getIndexHeader(char * passFile);
int getIndexCompatible(IndexHeader passHeader);

// Pack the given byte value into the BWT (pre-interleaved) at the specified index (4 per 32-bit word)
void packValue(bwt_t * passBWT, int64_t passSeqIdx, bwtint_t passValue);

// Unpack byte from from the BWT (pre-interleaved) at the specified index
ubyte_t unpackValue(bwt_t * passBWT, int64_t passSeqIdx);

// Obtain sequence length from the packed file
int64_t bwa_seq_len(const char *fn_pac);

// Populate the BWT from a packed file
bwt_t * bwt_pac2bwt(const char *fn_pac);

// 'pac2bwt' command entry point. (Note: bwt generated at this step CANNOT be used with BWA, bwtupdate required)
int command_pac2bwt(int argc, char *argv[]);

// Interleave occurrence counts into the BWT at specified interval for efficient search
void bwt_bwtupdate_core(bwt_t *bwt);

// 'bwtupdate' command entry point.
int command_bwtupdate(int argc, char *argv[]);

// 'bwt2sa' command entry point.
int command_bwt2sa(int argc, char *argv[]);

// 'index' command entry point.  Create protein file, pack, construct BWT, interleave, create SA, repack
int command_index(int argc, char *argv[]);

int64_t is_bwt(ubyte_t *T, int64_t n);

#endif /* BWTINDEX_H_ */
