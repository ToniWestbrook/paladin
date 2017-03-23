#ifndef PROTEIN_H_
#define PROTEIN_H_

#include "bwamem.h"
#include "bwtindex.h"

#define GFF_MAX_FIELD 8192

#define OUTPUT_TYPE_UNIPROT_SIMPLE 0
#define OUTPUT_TYPE_UNIPROT_FULL   1

#define ALIGN_FLAG_BRUTE_ORF  0x0001
#define ALIGN_FLAG_GEN_NT     0x0002
#define ALIGN_FLAG_KEEP_PRO   0x0004
#define ALIGN_FLAG_ADJUST_ORF 0x0008
#define ALIGN_FLAG_MANUAL_PRO 0x0010

typedef struct {
	unsigned long startIdx;
	unsigned long endIdx;
	short strand;
	short phase;
	short relFrame;
	char description[5000];
} CDS;

// Encoding
unsigned char encodeCodon(char * passSequence, int passStrand);
int convertToAA(char * passSequence, CDS * passCDS, int passTrans, char ** retSequence, unsigned long * retSize);

// ORF Detection
long getLastAlignedPos(long passLength, int passFrame);
long getLastAlignedOrfPos(long passLength, int passFrame, mem_opt_t * passOptions);
void addORFHistory(long * passHistory[2][6], long passHistorySize[6], unsigned long passIdx);
void compileORFHistory(long * passHistory[2][6], long passHistorySize[6], CDS * * retCDS, unsigned long * retCount);
int getSequenceORF(char * passSequence, unsigned long passLength, int passTrans, mem_opt_t * passOptions, CDS * * retCDS, unsigned long * retCount);

// Protein Creation
int getNextCDS(FILE * passFile, CDS * retCDS, unsigned long * retLine);
int writeIndexProtein(const char * passPrefix, const char * passProName, const char * passAnnName, IndexHeader passHeader);
int writeIndexCodingProtein(const char * passPrefix, const char * passProName, IndexHeader passHeader);
int writeIndexDirectProtein(const char * passPrefix, const char * passProName, IndexHeader passHeader);
int writeReadsProtein(const char * passPrefix, const char * passProName, mem_opt_t * passOptions);
int writeIndexTestProtein(const char * passPrefix, const char * proName);

#endif /* PROTEIN_H_ */
