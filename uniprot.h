// NOTE - Eventually separate this out into common reporting functionality and UniProt specific
// once we start supporting multiple report types/sources

#ifndef UNIPROT_H_
#define UNIPROT_H_

#include "bwamem.h"

#define UNIPROT_MAX_SUBMIT 5000
#define UNIPROT_MAX_ERROR 5
#define UNIPROT_BUFFER_GROW 50000000

#define UNIPROT_LIST_FULL 0
#define UNIPROT_LIST_GENES 1
#define UNIPROT_LIST_ORGANISM 2

#define UNIPROT_REFERENCE_SWISSPROT 1
#define UNIPROT_REFERENCE_UNIREF90 2

typedef struct {
	char * id;
	char * gene;
	char * organism;
	int numOccurrence;
	int totalQuality;
    int maxQuality;
} UniprotEntry;

typedef struct {
	UniprotEntry * entries;
	int entryCount;
	int unalignedCount;
} UniprotList;

typedef struct {
	char * buffer;
	int size;
	int capacity;
} CURLBuffer;

// Each pipeline run maintains its own list
extern UniprotList * uniprotPriEntryLists;
extern UniprotList * uniprotSecEntryLists;
extern int uniprotPriListCount;
extern int uniprotSecListCount;

// Rendering
void renderUniprotReport(int passType, int passPrimary, FILE * passStream, const char * passProxy);
void renderUniprotEntries(UniprotList * passList, int passType, FILE * passStream);
void renderNumberAligned(const mem_opt_t * passOptions);

// Population
int addUniprotList(worker_t * passWorker, int passSize, int passFull);
void cleanUniprotLists(UniprotList * passLists, int passPrimary);

// Support
UniprotList * getGlobalLists(int passPrimary);
int * getGlobalCount(int passPrimary);
void prepareUniprotReport(int passType, int passPrimary, UniprotList * passLists, CURLBuffer * passBuffer, const char * passProxy);
void prepareUniprotLists(UniprotList * retLists, int passPrimary);
void aggregateUniprotList(UniprotList * retList, int passListType, int passPrimary);
void joinOnlineLists(UniprotList * retList, char * passUniprotOutput);

// QSort Functions
int uniprotEntryCompareCommon (const void * passEntry1, const void * passEntry2);
int uniprotEntryCompareID (const void * passEntry1, const void * passEntry2);
int uniprotEntryCompareGene (const void * passEntry1, const void * passEntry2);
int uniprotEntryCompareOrganism (const void * passEntry1, const void * passEntry2);
int uniprotEntryCompareOnline (const void * passEntry1, const void * passEntry2);

// UniProt Interoperability
int cleanUniprotReference(int passReference, const char * passBase);
void cleanUniprotReferenceUniref(const char * passName, int passANN);
const char * downloadUniprotReference(int passReference, const char * passProxy);
void retrieveUniprotOnline(UniprotList * passList, CURLBuffer * retBuffer, const char * passProxy);
size_t receiveUniprotOutput(void * passString, size_t passSize, size_t passNum, void * retStream);
void initCURLBuffer(CURLBuffer * passBuffer, int passCapacity);
void resetCURLBuffer(CURLBuffer * passBuffer);
void freeCURLBuffer(CURLBuffer * passBuffer);

#endif /* UNIPROT_H_ */
