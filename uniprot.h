#ifndef UNIPROT_H_
#define UNIPROT_H_

#include "bwamem.h"

#define MAX_SUBMIT 300

#define UNIPROT_LIST_FULL 0
#define UNIPROT_LIST_ORGANISM 1
#define UNIPROT_LIST_GENES 2

typedef struct {
	char * * entries;
	int count;
} UniprotList;

static UniprotList * uniprotEntryLists = 0;
static int uniprotListCount = 0;

// Rendering
void renderUniprotReport(int passType);
void renderUniprotEntries(char * * passList, int * passCount, int passSize);

// Population
int addUniprotList(worker_t * passWorker, int passSize);
void cleanUniprotLists();

// Support
void prepareUniprotLists(char * * retLists[], int * retCounts[], int retSize[]);
int uniqueUniprotEntry(char * passValue, char * * passList, int * passCount, int passSize);
void sortUniprotEntries(char * * passList, int * passCount, int passSize);
void retrieveUniprotOnline(char * * passList, int * passCount, int passSize, char * * retOutput);
size_t receiveUniprotEntries(void * passString, size_t passSize, size_t passNum, void * retStream);
void joinOnlineLists();


#endif /* UNIPROT_H_ */
