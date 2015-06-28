/* Contact: Toni Westbrook <anthonyw@wildcats.unh.edu> */

#ifndef PROTEIN_H_
#define PROTEIN_H_

#include "bwamem.h"

#define GFF_MAX_FIELD 8192
#define INDEX_FLAG_NT 0x02
#define INDEX_FLAG_MF 0x04

extern unsigned char codon_aa_hash[64];

struct CDS {
	unsigned long startIdx;
	unsigned long endIdx;
	short strand;
	short phase;
	short relFrame;
	char description[5000];
};

// Encoding
void writeIndexHeader(FILE * passFilePtr, int passProtein, int passMulti);
char getIndexHeader(char * passFile);
unsigned char encodeCodon(unsigned char * passSequence, int passStrand);
int convertToAA(char * passSequence, struct CDS * passCDS, char ** retSequence, unsigned long * retSize);

// ORF Detection
long getLastORFPos(long passLength, int passFrame);
void addORFHistory(long * passHistory[2][6], long passHistorySize[6], unsigned long passIdx);
void compileORFHistory(long * passHistory[2][6], long passHistorySize[6], struct CDS * * retCDS, unsigned long * retCount);
int getSequenceORF(char * passSequence, unsigned long passLength, mem_opt_t * passOptions, char passIndexInfo, struct CDS * * retCDS, unsigned long * retCount);

// Protein Creation
int getNextCDS(FILE * passFile, struct CDS * retCDS, unsigned long * retLine);
int writeIndexProtein(const char * passPrefix, const char * passProName, const char * passAnnName);
int writeIndexMultiProtein(const char * passPrefix, const char * passProName, const char * passAnnName);
int writeIndexMultiCodingProtein(const char * passPrefix, const char * passProName);
int writeIndexDirectProtein(const char * passPrefix, const char * passProName);
int writeReadsProtein(const char * passPrefix, const char * passProName, mem_opt_t * passOptions, char passIndexInfo);

#endif /* PROTEIN_H_ */
