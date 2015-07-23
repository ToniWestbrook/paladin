#ifndef UNIPROT_H_
#define UNIPROT_H_

#include "bwamem.h"

#define UNIPROT_MAX_SUBMIT 10000
#define UNIPROT_MAX_ERROR 5
#define UNIPROT_BUFFER_GROW 50000000

#define UNIPROT_LIST_FULL 0
#define UNIPROT_LIST_GENES 1
#define UNIPROT_LIST_ORGANISM 2

typedef struct {
	char * id;
	char * gene;
	char * organism;
	int numOccurrence;
} UniprotEntry;

typedef struct {
	UniprotEntry * entries;
	int count;
} UniprotList;

typedef struct {
	char * buffer;
	int size;
	int capacity;
} CURLBuffer;

// Each pipeline run maintains its own list
extern UniprotList * uniprotEntryLists;
extern int uniprotListCount;

// Rendering
void renderUniprotReport(int passType);
void renderUniprotEntries(UniprotList * passList, int passType);

// Population
int addUniprotList(worker_t * passWorker, int passSize);
void cleanUniprotLists(UniprotList * passLists);

// Support
void prepareUniprotLists(UniprotList * retLists);
void aggregateUniprotList(UniprotList * retList, int passListType);
void joinOnlineLists(UniprotList * retList, char * passUniprotOutput);

// QSort Functions
int uniprotEntryCompareID (const void * passEntry1, const void * passEntry2);
int uniprotEntryCompareGene (const void * passEntry1, const void * passEntry2);
int uniprotEntryCompareOrganism (const void * passEntry1, const void * passEntry2);
int uniprotEntryCompareOnline (const void * passEntry1, const void * passEntry2);

// UniProt Interoperability
const char * downloadUniprotReference(int passReference);
void retrieveUniprotOnline(UniprotList * passList, CURLBuffer * retBuffer);
size_t receiveUniprotOutput(void * passString, size_t passSize, size_t passNum, void * retStream);
void initCURLBuffer(CURLBuffer * passBuffer, int passCapacity);
void resetCURLBuffer(CURLBuffer * passBuffer);
void freeCURLBuffer(CURLBuffer * passBuffer);


#endif /* UNIPROT_H_ */
