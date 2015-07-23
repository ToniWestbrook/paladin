#include <stdlib.h>
#include <string.h>
#include <curl/curl.h>
#include "uniprot.h"
#include "protein.h"

UniprotList * uniprotEntryLists = 0;
int uniprotListCount = 0;

const char * downloadNames[] = {"uniprot_sprot.fasta.gz", 
			  "uniprot_trembl.fasta.gz"};
const char * downloadURLs[] = {"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
			 "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"};

void renderUniprotReport(int passType) {
	UniprotList uniprotLists[3];
	CURLBuffer tempBuffer;

	// Aggregate and sort lists
	prepareUniprotLists(uniprotLists);

	if (uniprotLists[0].count == 0) {
		printf("No entries to report\n");
		return;
	}

	// Render requested report
	switch (passType) {
		case OUTPUT_TYPE_UNIPROT_SIMPLE:
			// Sort all local lists
			qsort(uniprotLists[UNIPROT_LIST_FULL].entries, uniprotLists[UNIPROT_LIST_FULL].count, sizeof(UniprotEntry), uniprotEntryCompareID);
			qsort(uniprotLists[UNIPROT_LIST_GENES].entries, uniprotLists[UNIPROT_LIST_GENES].count, sizeof(UniprotEntry), uniprotEntryCompareGene);
			qsort(uniprotLists[UNIPROT_LIST_ORGANISM].entries, uniprotLists[UNIPROT_LIST_ORGANISM].count, sizeof(UniprotEntry), uniprotEntryCompareOrganism);

			// Render output
			printf("Count\tUniProtKB\n");
			renderUniprotEntries(uniprotLists + UNIPROT_LIST_FULL, UNIPROT_LIST_FULL);
			printf("\n\nCount\tGene\n");
			renderUniprotEntries(uniprotLists + UNIPROT_LIST_GENES, UNIPROT_LIST_GENES);
			printf("\n\nCount\tOrganism\n");
			renderUniprotEntries(uniprotLists + UNIPROT_LIST_ORGANISM, UNIPROT_LIST_ORGANISM);

			break;

		case OUTPUT_TYPE_UNIPROT_FULL:
			// Submit entries to UniProt and retrieve full information
			retrieveUniprotOnline(uniprotLists + UNIPROT_LIST_FULL, &tempBuffer);
			joinOnlineLists(uniprotLists + UNIPROT_LIST_FULL, tempBuffer.buffer);
			qsort(uniprotLists[UNIPROT_LIST_FULL].entries, uniprotLists[UNIPROT_LIST_FULL].count, sizeof(UniprotEntry), uniprotEntryCompareID);

			// Render output
			printf("Count\tUniProtKB\tID\tOrganism\tProtein Names\tGenes\tPathway\tFeatures\tGene Ontology\tReviewd\tExistence\tComments\n");
			renderUniprotEntries(uniprotLists + UNIPROT_LIST_FULL, UNIPROT_LIST_FULL);
			freeCURLBuffer(&tempBuffer);
			break;
	}

	cleanUniprotLists(uniprotLists);
}

const char * downloadUniprotReference(int passReference) {
	CURL * curlHandle;
	CURLcode curlResult;
	FILE * fileHandle;
	const char * retFile;

	curlHandle = curl_easy_init();
	fileHandle = fopen(downloadNames[passReference], "w");
	retFile = downloadNames[passReference];

	curl_easy_setopt(curlHandle, CURLOPT_URL, downloadURLs[passReference]);
	curl_easy_setopt(curlHandle, CURLOPT_WRITEDATA, fileHandle) ;

	fprintf(stderr, "[M::%s] Downloading %s...\n", __func__, downloadURLs[passReference]);

	curlResult = curl_easy_perform(curlHandle);

	if (curlResult != CURLE_OK) {
		fprintf(stderr, "ERROR: %s\n", curl_easy_strerror(curlResult));
		unlink(downloadNames[passReference]);
		retFile = "";
	}
		
	curl_easy_cleanup(curlHandle);
	fclose(fileHandle);

	return retFile;
}

void retrieveUniprotOnline(UniprotList * passList, CURLBuffer * retBuffer) {
	int entryIdx, queryIdx, parseIdx, errorIdx, queryCount;
	CURL * curlHandle;
	CURLcode curlResult;
	CURLBuffer tempBuffer;
	char queryString[UNIPROT_MAX_SUBMIT * 50], jobID[50];
	char * httpString;

	// Init structures
	curlHandle = curl_easy_init();
	initCURLBuffer(retBuffer, UNIPROT_BUFFER_GROW);
	initCURLBuffer(&tempBuffer, UNIPROT_BUFFER_GROW);
	queryCount = (passList->count < UNIPROT_MAX_SUBMIT) ? passList->count : UNIPROT_MAX_SUBMIT;

	for (entryIdx = 0 ; entryIdx < passList->count ; ) {
		// Build query string
		queryString[0] = 0;
		for (queryIdx = 0 ; (queryIdx < queryCount) && (entryIdx < passList->count) ; entryIdx++) {
			for (parseIdx = 0 ; parseIdx < strlen(passList->entries[entryIdx].id) ; parseIdx++) {
				if (passList->entries[entryIdx].id[parseIdx] == '_') {
					sprintf(queryString, "%s%s ", queryString, passList->entries[entryIdx].id);
					queryIdx++;
					break;
				}
			}
		}

		httpString = curl_easy_escape(curlHandle, queryString, 0);
		sprintf(queryString, "uploadQuery=%s&format=job&from=ACC+ID&to=ACC&landingPage=false", httpString);
		curl_free(httpString);

		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Submitted %d of %d entries to UniProt...\n", __func__, entryIdx, passList->count);
		}

		// Restart a limited number of times if errors encountered
		for (errorIdx = 0 ; errorIdx < UNIPROT_MAX_ERROR ; errorIdx++) {
			// Stage 1 - Submit query for processing
			curl_easy_setopt(curlHandle, CURLOPT_URL, "http://www.uniprot.org/uploadlists/");
			curl_easy_setopt(curlHandle, CURLOPT_POSTFIELDS, queryString);
			curl_easy_setopt(curlHandle, CURLOPT_FOLLOWLOCATION, 1L);
			curl_easy_setopt(curlHandle, CURLOPT_WRITEFUNCTION, receiveUniprotOutput);
			curl_easy_setopt(curlHandle, CURLOPT_WRITEDATA, &tempBuffer);

			resetCURLBuffer(&tempBuffer);
			curlResult = curl_easy_perform(curlHandle);

			if (curlResult != CURLE_OK) {
				fprintf(stderr, "ERROR: %s\n", curl_easy_strerror(curlResult));
				continue;
			}

			// Stage 2 - Wait for results
			sprintf(jobID, "%s", tempBuffer.buffer);
			sprintf(queryString, "http://www.uniprot.org/jobs/%s.stat", jobID);

			resetCURLBuffer(&tempBuffer);
			curl_easy_setopt(curlHandle, CURLOPT_URL, queryString);

			while (strcmp(tempBuffer.buffer, "COMPLETED") != 0) {
				sleep(1);
				resetCURLBuffer(&tempBuffer);
				curlResult = curl_easy_perform(curlHandle);
				if (curlResult != CURLE_OK) break;
			}

			if (curlResult != CURLE_OK) {
				fprintf(stderr, "ERROR: %s\n", curl_easy_strerror(curlResult));
				continue;
			}

			// Stage 3 - Retrieve Results
			sprintf(queryString, "query=job:%s&format=tab&columns=entry%%20name,id,organism,protein%%20names,genes,pathway,features,go,reviewed,existence,comments", jobID);
			curl_easy_setopt(curlHandle, CURLOPT_URL, "http://www.uniprot.org/uniprot/");
			curl_easy_setopt(curlHandle, CURLOPT_WRITEDATA, retBuffer);

			curlResult = curl_easy_perform(curlHandle);

			if (curlResult != CURLE_OK) {
				fprintf(stderr, "ERROR: %s\n", curl_easy_strerror(curlResult));
				continue;
			}

			break;
		}
	}

	curl_easy_cleanup(curlHandle);
	freeCURLBuffer(&tempBuffer);
}

void renderUniprotEntries(UniprotList * passList, int passType) {
	int entryIdx;


	for (entryIdx = 0 ; entryIdx < passList->count ; entryIdx++) {
		switch(passType) {
			case UNIPROT_LIST_FULL:
				printf("%d\t%s\n", passList->entries[entryIdx].numOccurrence, passList->entries[entryIdx].id);
				break;
			case UNIPROT_LIST_GENES:
				printf("%d\t%s\n", passList->entries[entryIdx].numOccurrence, passList->entries[entryIdx].gene);
				break;
			case UNIPROT_LIST_ORGANISM:
				printf("%d\t%s\n", passList->entries[entryIdx].numOccurrence, passList->entries[entryIdx].organism);
				break;
		}
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
	uniprotEntryLists = realloc(uniprotEntryLists, (uniprotListCount + 1) * sizeof(UniprotList));
	uniprotEntryLists[uniprotListCount].entries = calloc(listSize, sizeof(UniprotEntry));
	uniprotEntryLists[uniprotListCount].count = listSize;

	// Populate list
	for (addIdx = 0, entryIdx = 0 ; entryIdx < passSize ; entryIdx++) {
		for (alnIdx = 0 ; alnIdx < passWorker->regs[entryIdx].n ; alnIdx++) {

			// Extract ID and entry
			refID = passWorker->regs[entryIdx].a[alnIdx].rid;
			uniprotEntry = passWorker->bns->anns[refID].name;

			// Strip sequence and frame info for nucleotide references
			if (passWorker->opt->indexFlag & INDEX_FLAG_NT)
			for (parseIdx = 2 ; (parseIdx > 0) && (*uniprotEntry != 0) ; uniprotEntry++) {
				if (*uniprotEntry == ':') parseIdx--;
			}

			// Strip initial IDs
			for (parseIdx = 2 ; (parseIdx > 0) && (*uniprotEntry != 0) ; uniprotEntry++) {
				if (*uniprotEntry == '|') parseIdx--;
			}

			// Strip description
			for (parseIdx = 0 ; uniprotEntry[parseIdx] != 0 ; parseIdx++) {
				if (uniprotEntry[parseIdx] == ' ') {
					uniprotEntry[parseIdx] = 0;
					break;
				}
			}

			// Full ID
			uniprotEntryLists[uniprotListCount].entries[addIdx].id = malloc(strlen(uniprotEntry) + 1);
			uniprotEntryLists[uniprotListCount].entries[addIdx].numOccurrence = 1;
			sprintf(uniprotEntryLists[uniprotListCount].entries[addIdx].id, "%s", uniprotEntry);

			// Gene/organism
			for (parseIdx = 0 ; parseIdx < strlen(uniprotEntry) ; parseIdx++) {
				if (*(uniprotEntry + parseIdx) == '_') {
					uniprotEntryLists[uniprotListCount].entries[addIdx].gene = malloc(parseIdx + 1);
					sprintf(uniprotEntryLists[uniprotListCount].entries[addIdx].gene, "%.*s", parseIdx, uniprotEntry);
					uniprotEntryLists[uniprotListCount].entries[addIdx].organism = malloc(strlen(uniprotEntry + parseIdx) + 1);
					sprintf(uniprotEntryLists[uniprotListCount].entries[addIdx].organism, "%s", uniprotEntry + parseIdx + 1);
					break;
				}
			}

			addIdx++;
		}
	}

	return uniprotListCount++;
}

void cleanUniprotLists(UniprotList * passLists) {
	int listIdx, entryIdx;

	// Global list 0 contains pointers to all allocated strings, so delete from there. Do not repeat with other lists
	for (entryIdx = 0 ; entryIdx < uniprotEntryLists[0].count ; entryIdx++) {
		free(uniprotEntryLists[0].entries[entryIdx].id);
		free(uniprotEntryLists[0].entries[entryIdx].gene);
		free(uniprotEntryLists[0].entries[entryIdx].organism);
	}

	free(uniprotEntryLists[0].entries);
	free(uniprotEntryLists);

	// Free local lists
	for (listIdx = 0 ; listIdx < 3 ; listIdx++) {
		free(passLists[listIdx].entries);
	}
}

size_t receiveUniprotOutput(void * passString, size_t passSize, size_t passNum, void * retStream) {
	CURLBuffer * currentBuffer;

	currentBuffer = (CURLBuffer *) retStream;

	// Grow receive buffer if addition is greater than capacity
	if (currentBuffer->size + (int)(passSize * passNum) >= currentBuffer->capacity) {
		currentBuffer->buffer = realloc(currentBuffer->buffer, currentBuffer->capacity + UNIPROT_BUFFER_GROW);
		currentBuffer->capacity += UNIPROT_BUFFER_GROW;
        }

	// Concatenate results
	memcpy(currentBuffer->buffer + currentBuffer->size, passString, (int) (passSize * passNum));
	currentBuffer->size += (int) (passSize * passNum);
	currentBuffer->buffer[currentBuffer->size] = 0;

	return (size_t)(passSize * passNum);
}


void prepareUniprotLists(UniprotList * retLists) {
	int listIdx, entryIdx, localIdx;
	int maxEntries, totalSize;

	// Calculate maximum number of UniProt entries
	for (maxEntries = 0, listIdx = 0 ; listIdx < uniprotListCount ; listIdx++) {
		maxEntries += uniprotEntryLists[listIdx].count;
	}

	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Aggregating %d entries for UniProt report\n", __func__, maxEntries);
	}

	if (maxEntries == 0) return;

	// Join each pipeline's individual list into one master list (at first list)
	for (listIdx = 0, totalSize = 0 ; listIdx < uniprotListCount ; listIdx++) {
		totalSize += uniprotEntryLists[listIdx].count;
	}

	uniprotEntryLists[0].entries = realloc(uniprotEntryLists[0].entries, sizeof(UniprotEntry) * totalSize);

	for (listIdx = 1, entryIdx = uniprotEntryLists[0].count ; listIdx < uniprotListCount ; listIdx++) {
		memcpy(uniprotEntryLists[0].entries + entryIdx, uniprotEntryLists[listIdx].entries, sizeof(UniprotEntry) * uniprotEntryLists[listIdx].count);
		entryIdx += uniprotEntryLists[listIdx].count;
		free(uniprotEntryLists[listIdx].entries);
		uniprotEntryLists[listIdx].count = 0;
	}

	uniprotEntryLists[0].count = totalSize;

	// Aggregate each local list
	for (localIdx = 0 ; localIdx < 3 ; localIdx++) {
		aggregateUniprotList(retLists + localIdx, localIdx);
	}
}

void joinOnlineLists(UniprotList * retList, char * passUniprotOutput) {
	char * * lineIndices;
	int parseIdx, entryIdx, lineIdx;
	int lineCount, outputSize, matchValue;

	// Count number of lines in output
	outputSize = strlen(passUniprotOutput);

	for (lineCount = 0, parseIdx = 0 ; parseIdx < outputSize ; parseIdx++) {
		if (passUniprotOutput[parseIdx] == '\n') lineCount++;
	}
	if (passUniprotOutput[parseIdx - 1] != '\n') lineCount++;

	// Index each line
	lineIndices = malloc(lineCount * sizeof(char *));
	lineIndices[0] = passUniprotOutput;

	for (lineIdx = 1, parseIdx = 0 ; parseIdx < outputSize ; parseIdx++) {
		// If EOL found, change to NULL and check for next valid line
		if (passUniprotOutput[parseIdx] == '\n') {
			passUniprotOutput[parseIdx] = 0;
			if (parseIdx < outputSize - 1) {
				lineIndices[lineIdx++] = passUniprotOutput + parseIdx + 1;
			}
		}
	}

	// Sort (hopefully temporarily)
	qsort(lineIndices, lineCount, sizeof(char *), uniprotEntryCompareOnline);

	// Now cross reference/join - both are ordered, so skip in lexicographical order if match not found
	for (entryIdx = 0, lineIdx = 1 ; (entryIdx < retList->count) && (lineIdx < lineCount)  ; ) {
		matchValue = strncmp(retList->entries[entryIdx].id, lineIndices[lineIdx], strlen(retList->entries[entryIdx].id));

		if (matchValue == 0) {
			// Do not free existing string - global list handles this during cleanup
			retList->entries[entryIdx++].id = lineIndices[lineIdx++];
		}
		else if (matchValue > 0) lineIdx++;
		else entryIdx++; 
	}

	free(lineIndices);
}

void aggregateUniprotList(UniprotList * retList, int passListType) {
	int entryIdx;
	int memberOffset;
	char * srcBase, * dstBase;

	// First sort full list
	switch (passListType) {
		case UNIPROT_LIST_FULL:
			qsort(uniprotEntryLists[0].entries, uniprotEntryLists[0].count, sizeof(UniprotEntry), uniprotEntryCompareID);
			memberOffset = offsetof(UniprotEntry, id);
			break;
		case UNIPROT_LIST_GENES:
			qsort(uniprotEntryLists[0].entries, uniprotEntryLists[0].count, sizeof(UniprotEntry), uniprotEntryCompareGene);
			memberOffset = offsetof(UniprotEntry, gene);
			break;
		case UNIPROT_LIST_ORGANISM:
			qsort(uniprotEntryLists[0].entries, uniprotEntryLists[0].count, sizeof(UniprotEntry), uniprotEntryCompareOrganism);
			memberOffset = offsetof(UniprotEntry, organism);
			break;
		default: memberOffset = 0; break;
	}

	// Count entry occurrences and aggregate into return list
	retList->entries = calloc(uniprotEntryLists[0].count, sizeof(UniprotEntry));
	retList->entries[0].id = uniprotEntryLists[0].entries[0].id;
	retList->entries[0].gene = uniprotEntryLists[0].entries[0].gene;
	retList->entries[0].organism = uniprotEntryLists[0].entries[0].organism;
	retList->count = 0;

	for (entryIdx = 0 ; entryIdx < uniprotEntryLists[0].count ; entryIdx++) {
		srcBase = (char *) (uniprotEntryLists[0].entries + entryIdx);
		dstBase = (char *) (retList->entries + retList->count);

		// If current global entry doesn't match local entry, create new local entry and copy pointer to relevant string info
		if (strcmp(*((char * *) (srcBase + memberOffset)), *((char * *) (dstBase + memberOffset))) != 0) {
			dstBase = (char *) (retList->entries + ++retList->count);
			*((char * *) (dstBase + memberOffset)) = *((char * *) (srcBase + memberOffset));
		}

		retList->entries[retList->count].numOccurrence++;
	}

	retList->count++;
}

int uniprotEntryCompareID (const void * passEntry1, const void * passEntry2) {
	if (((UniprotEntry *)passEntry1)->numOccurrence == ((UniprotEntry *)passEntry2)->numOccurrence) {
		return (strcmp(((UniprotEntry *)passEntry1)->id, ((UniprotEntry *)passEntry2)->id));
	}
	else {
		return ((UniprotEntry *)passEntry2)->numOccurrence - ((UniprotEntry *)passEntry1)->numOccurrence;
	}
}

int uniprotEntryCompareGene (const void * passEntry1, const void * passEntry2) {
	if (((UniprotEntry *)passEntry1)->numOccurrence == ((UniprotEntry *)passEntry2)->numOccurrence) {
		return (strcmp(((UniprotEntry *)passEntry1)->gene, ((UniprotEntry *)passEntry2)->gene));
	}
	else {
		return ((UniprotEntry *)passEntry2)->numOccurrence - ((UniprotEntry *)passEntry1)->numOccurrence;
	}
}

int uniprotEntryCompareOrganism (const void * passEntry1, const void * passEntry2) {
	if (((UniprotEntry *)passEntry1)->numOccurrence == ((UniprotEntry *)passEntry2)->numOccurrence) {
		return (strcmp(((UniprotEntry *)passEntry1)->organism, ((UniprotEntry *)passEntry2)->organism));
	}
	else {
		return ((UniprotEntry *)passEntry2)->numOccurrence - ((UniprotEntry *)passEntry1)->numOccurrence;
	}
}

int uniprotEntryCompareOnline (const void * passEntry1, const void * passEntry2) {
	return (strcmp(*((char * *)passEntry1), *((char * *)passEntry2)));
}

void initCURLBuffer(CURLBuffer * passBuffer, int passCapacity) {
	passBuffer->buffer = calloc(1, passCapacity);
	passBuffer->capacity = passCapacity;
	passBuffer->size = 0;
	passBuffer->buffer[0] = 0;
}

void resetCURLBuffer(CURLBuffer * passBuffer) {
	passBuffer->size = 0;
	passBuffer->buffer[0] = 0;
}

void freeCURLBuffer(CURLBuffer * passBuffer) {
	free(passBuffer->buffer);
}
