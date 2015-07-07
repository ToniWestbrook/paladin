#include <stdlib.h>
#include <string.h>
#include <curl/curl.h>
#include "uniprot.h"
#include "protein.h"

void renderUniprotReport(int passType) {
	char * * uniprotLists[3];
	int * uniprotCounts[3];
	int uniprotSize[3];
	char * tempOutput;
	int entriesRetrieved;

	// Aggregate and sort lists
	prepareUniprotLists(uniprotLists, uniprotCounts, uniprotSize);

	// Render requested report
	switch (passType) {
		case OUTPUT_TYPE_UNIPROT_SIMPLE:
			printf("Count\tUniProtKB\n");
			renderUniprotEntries(uniprotLists[UNIPROT_LIST_FULL], uniprotCounts[UNIPROT_LIST_FULL], uniprotSize[UNIPROT_LIST_FULL]);

			printf("\n\nCount\tGene\n");
			renderUniprotEntries(uniprotLists[UNIPROT_LIST_GENES], uniprotCounts[UNIPROT_LIST_GENES], uniprotSize[UNIPROT_LIST_GENES]);

			printf("\n\nCount\tOrganism\n");
			renderUniprotEntries(uniprotLists[UNIPROT_LIST_ORGANISM], uniprotCounts[UNIPROT_LIST_ORGANISM], uniprotSize[UNIPROT_LIST_ORGANISM]);

			break;

		case OUTPUT_TYPE_UNIPROT_FULL:
			printf("Count\tUniProtKB\tID\tOrganism\tProtein Names\tGenes\tPathway\tFeatures\tGene Ontology\tReviewd\tExistence\tComments\n");
			retrieveUniprotOnline(uniprotLists[UNIPROT_LIST_FULL], uniprotCounts[UNIPROT_LIST_FULL], uniprotSize[UNIPROT_LIST_FULL], &tempOutput);
			joinOnlineLists(uniprotLists[UNIPROT_LIST_FULL], uniprotSize[UNIPROT_LIST_FULL], tempOutput);
			renderUniprotEntries(uniprotLists[UNIPROT_LIST_FULL], uniprotCounts[UNIPROT_LIST_FULL], uniprotSize[UNIPROT_LIST_FULL]);
			break;
	}
}

void retrieveUniprotOnline(char * * passList, int * passCount, int passSize, char * * retOutput) {
	int listIdx, printIdx, parseIdx, queryCount;
	CURL * curlHandle;
	CURLcode curlResult;
	char queryString[passSize * 50];
	char * httpString;

	listIdx = 0;
	*retOutput = calloc(1, 1);
	curlHandle = curl_easy_init();

	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Submitting %d entries to UniProt (this may take a while)...", __func__, passSize);
	}

	// Build query string
	queryCount = (passSize < MAX_SUBMIT) ? passSize : MAX_SUBMIT;

	while (listIdx < passSize) {
		queryString[0] = 0;

		for (printIdx = 0 ; (printIdx < queryCount) && (listIdx < passSize) ; listIdx++) {
			for (parseIdx = 0 ; parseIdx < strlen(passList[listIdx]) ; parseIdx++) {
				if (passList[listIdx][parseIdx] == '_') {
					sprintf(queryString, "%s%s or ", queryString, passList[listIdx]);
					printIdx++;
					break;
				}
			}
		}

		// Finalize query string
		if (queryString[0] != 0) {
			queryString[strlen(queryString) - 4] = 0;
		}

		httpString = curl_easy_escape(curlHandle, queryString, 0);
		sprintf(queryString, "query=%s&format=tab&columns=entry%%20name,id,organism,protein%%20names,genes,pathway,features,go,reviewed,existence,comments", httpString);
		curl_free(httpString);

		curl_easy_setopt(curlHandle, CURLOPT_URL, "http://www.uniprot.org/uniprot/");
		curl_easy_setopt(curlHandle, CURLOPT_POSTFIELDS, queryString);
		curl_easy_setopt(curlHandle, CURLOPT_FOLLOWLOCATION, 1L);

		curl_easy_setopt(curlHandle, CURLOPT_WRITEFUNCTION, receiveUniprotEntries);
		curl_easy_setopt(curlHandle, CURLOPT_WRITEDATA, retOutput);
		curlResult = curl_easy_perform(curlHandle);

		if (curlResult != CURLE_OK) {
			fprintf(stderr, "ERROR: %s\n", curl_easy_strerror(curlResult));
		}

		if (bwa_verbose >= 3) {
			fprintf(stderr, ".");
		}
	}

	curl_easy_cleanup(curlHandle);

	if (bwa_verbose >= 3) {
		fprintf(stderr, "\n");
	}
}

void renderUniprotEntries(char * * passList, int * passCount, int passSize) {
	int listIdx;

	for (listIdx = 0 ; listIdx < passSize ; listIdx++) {
		printf("%d\t%s\n", passCount[listIdx], passList[listIdx]);
	}
}

int addUniprotList(worker_t * passWorker, int passSize) {
	int entryIdx, alnIdx, addIdx, parseIdx;
	int refID, listSize;
	char * uniprotEntry;

	// Calculate total list size
	for (listSize = 0 , entryIdx = 0 ; entryIdx < passSize ; entryIdx++) {
		listSize += passWorker->regs[entryIdx].n;
	}

	// Create list
	uniprotEntryLists = realloc(uniprotEntryLists, uniprotListCount);
	uniprotEntryLists[uniprotListCount].entries = calloc(listSize, sizeof(char *));
	uniprotEntryLists[uniprotListCount].count = listSize;

	// Populate list
	for (addIdx = 0, entryIdx = 0 ; entryIdx < passSize ; entryIdx++) {
		for (alnIdx = 0 ; alnIdx < passWorker->regs[entryIdx].n ; alnIdx++) {

			// Extract ID and entry
			refID = passWorker->regs[entryIdx].a[alnIdx].rid;
			uniprotEntry = passWorker->bns->anns[refID].name;

			// Strip sequence and frame info
			for (parseIdx = 2 ; (parseIdx > 0) && (*uniprotEntry != 0) ; uniprotEntry++) {
				if (*uniprotEntry == ':') parseIdx--;
			}

			uniprotEntryLists[uniprotListCount].entries[addIdx] = malloc(strlen(uniprotEntry) + 1);
			sprintf(uniprotEntryLists[uniprotListCount].entries[addIdx], "%s", uniprotEntry);
			addIdx++;
		}

		free(passWorker->regs[entryIdx].a);
	}

	return uniprotListCount++;
}

void cleanUniprotLists() {
	int listIdx, entryIdx;

	for (listIdx = 0 ; listIdx < uniprotListCount ; listIdx++) {
		for (entryIdx = 0 ; entryIdx < uniprotEntryLists[listIdx].count ; entryIdx++) {
			free(uniprotEntryLists[listIdx].entries[entryIdx]);
		}
		free(uniprotEntryLists[listIdx].entries);
	}

	free(uniprotEntryLists);
}

int uniqueUniprotEntry(char * passValue, char * * passList, int * passCount, int passSize) {
	int entryIdx, retVal;

	retVal = 0;

	for (entryIdx = 0 ; entryIdx < passSize ; entryIdx++) {
		// If we've hit the end of the list, add entry
		if (passList[entryIdx] == 0) {
			passList[entryIdx] = malloc(strlen(passValue) + 1);
			sprintf(passList[entryIdx], "%s", passValue);
			retVal++;
			break;
		}

		// Check for match
		if (strcmp(passValue, passList[entryIdx]) == 0) break;
	}

	// Increase count
	passCount[entryIdx]++;

	return retVal;
}

void sortUniprotEntries(char * * passList, int * passCount, int passSize) {
	int idx1, idx2;
	int tempCount;
	char * tempValue;

	// Sort values
	for (idx1 = passSize - 1 ; idx1 > 0 ; idx1--) {
		for (idx2 = 0 ; idx2 < idx1 ; idx2++) {
			// First sort descending by size, then ascending by name
			if ((passCount[idx2] < passCount[idx2 + 1]) ||
				((passCount[idx2] == passCount[idx2 + 1]) && (strcmp(passList[idx2], passList[idx2 + 1]) > 0))) {
				tempCount = passCount[idx2 + 1];
				tempValue = passList[idx2 + 1];
				passCount[idx2 + 1] = passCount[idx2];
				passList[idx2 + 1] = passList[idx2];
				passCount[idx2] = tempCount;
				passList[idx2] = tempValue;
			}
		}
	}
}

size_t receiveUniprotEntries(void * passString, size_t passSize, size_t passNum, void * retStream) {
	char * tempPtr;

	// Concatenate previous results to new
	tempPtr = malloc(strlen(*((char * *) retStream)) + (passSize * passNum) + 1);
	sprintf(tempPtr, "%s%.*s", *((char * *) retStream), (int) (passSize * passNum), (char *) passString);

	// Return results
	free(*((char * *) retStream));
	*((char * *)retStream) = tempPtr;

	return (size_t)(passSize * passNum);
}


void prepareUniprotLists(char * * retLists[], int * retCounts[], int retSize[]) {
	int listIdx, entryIdx, parseIdx, localIdx;
	int maxEntries;
	char * uniprotEntry;

	// Calculate maximum number of UniProt entries
	for (maxEntries = 0, listIdx = 0 ; listIdx < uniprotListCount ; listIdx++) {
		maxEntries += uniprotEntryLists[listIdx].count;
	}

	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Aggregating %d entries for UniProt report\n", __func__, maxEntries);
	}

	// Initialize local
	for (localIdx = 0 ; localIdx < 3 ; localIdx++) {
		retSize[localIdx] = 0;
		retLists[localIdx] = calloc(maxEntries, sizeof(char *));
		retCounts[localIdx] = calloc(maxEntries, sizeof(int));
	}

	// Add unique entries to lists
	for (listIdx = 0 ; listIdx < uniprotListCount ; listIdx++) {
		for (entryIdx = 0 ; entryIdx < uniprotEntryLists[listIdx].count ; entryIdx++) {
			uniprotEntry = uniprotEntryLists[listIdx].entries[entryIdx];

			// Add unique entry to list
			retSize[UNIPROT_LIST_FULL] += uniqueUniprotEntry(uniprotEntry, retLists[UNIPROT_LIST_FULL], retCounts[UNIPROT_LIST_FULL], maxEntries);

			// Parse gene/organism
			for (parseIdx = 0 ; parseIdx < strlen(uniprotEntry) ; parseIdx++) {
				if (*(uniprotEntry + parseIdx) == '_') {
					*(uniprotEntry + parseIdx) = 0;
					break;
				}
			}

			retSize[UNIPROT_LIST_GENES] += uniqueUniprotEntry(uniprotEntry, retLists[UNIPROT_LIST_GENES], retCounts[UNIPROT_LIST_GENES], maxEntries);
			retSize[UNIPROT_LIST_ORGANISM] += uniqueUniprotEntry(uniprotEntry + parseIdx + 1, retLists[UNIPROT_LIST_ORGANISM], retCounts[UNIPROT_LIST_ORGANISM], maxEntries);
		}
	}

	// Sort lists
	for (localIdx = 0 ; localIdx < 3 ; localIdx++) {
		sortUniprotEntries(retLists[localIdx], retCounts[localIdx], retSize[localIdx]);
	}
}

void joinOnlineLists(char * * passFullList, int passListSize, char * passUniprotOutput) {
	char * * lineIndices;
	int parseIdx, entryIdx, listIdx, entryCount, outputSize;

	// Count number of lines in output
	outputSize = strlen(passUniprotOutput);

	for (entryCount = 0, parseIdx = 0 ; parseIdx < outputSize ; parseIdx++) {
		if (passUniprotOutput[parseIdx] == '\n') entryCount++;
	}
	if (passUniprotOutput[parseIdx - 1] != '\n') entryCount++;

	// Index each line
	lineIndices = malloc(entryCount * sizeof(char *));
	lineIndices[0] = passUniprotOutput;

	for (entryIdx = 1, parseIdx = 0 ; parseIdx < outputSize ; parseIdx++) {
		// If EOL found, change to NULL and check for next valid line
		if (passUniprotOutput[parseIdx] == '\n') {
			passUniprotOutput[parseIdx] = 0;
			if (parseIdx < outputSize - 1) {
				lineIndices[entryIdx++] = passUniprotOutput + parseIdx + 1;
			}
		}
	}

	// Now cross reference lists and join into original
	for (listIdx = 0 ; listIdx < passListSize ; listIdx++) {
		for (entryIdx = 0 ; entryIdx < entryCount ; entryIdx++) {
			if (strncmp(passFullList[listIdx], lineIndices[entryIdx], strlen(passFullList[listIdx])) == 0) {
				free(passFullList[listIdx]);
				passFullList[listIdx] = lineIndices[entryIdx];
				break;
			}
		}
	}
}
