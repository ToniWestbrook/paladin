/* Contact: Toni Westbrook <anthonyw@wildcats.unh.edu> */

#ifndef PROTEIN_H_
#define PROTEIN_H_

#define GFF_MAX_FIELD 8192

extern unsigned char codon_aa_hash[64];

struct CDS {
	unsigned long startIdx;
	unsigned long endIdx;
	short strand;
	short phase;
	char description[4096];
};

// Encoding
unsigned char encodeCodon(unsigned char * passSequence, int passStrand);
int convertToAA(char * passSequence, struct CDS * passCDS, char ** retSequence, unsigned long * retSize);

// ORF Detection
long getLastORFPos(long passLength, int passFrame);
void addORFHistory(long * passHistory[2][6], long passHistorySize[6], unsigned long passIdx);
void compileORFHistory(long * passHistory[2][6], long passHistorySize[6], struct CDS * * retCDS, unsigned long * retCount);
int getSequenceORF(char * passSequence, unsigned long passLength, struct CDS * * retCDS, unsigned long * retCount);

// Protein Creation
int getNextCDS(FILE * passFile, struct CDS * retCDS, unsigned long * retLine);
int writeIndexProtein(const char * passPrefix, const char * passProName, const char * passAnnName);
int writeIndexMultiProtein(const char * passPrefix, const char * passProName, const char * passAnnName);
int writeReadsProtein(const char * passPrefix, const char * passProName);

#endif /* PROTEIN_H_ */
