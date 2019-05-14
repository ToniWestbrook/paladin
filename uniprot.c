#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <curl/curl.h>
#include <unistd.h>
#include <zlib.h>
#include "uniprot.h"
#include "main.h"
#include "protein.h"
#include "utils.h"

UniprotList * uniprotPriEntryLists = 0;
UniprotList * uniprotSecEntryLists = 0;
int uniprotPriListCount = 0;
int uniprotSecListCount = 0;

// For code clarity, 0 position reserved for non-Uniprot reference
const char * downloadNames[] = {"",
		"uniprot_sprot.fasta.gz",
		"uniref90.fasta.gz"};
const char * downloadURLs[] = {"",
		"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
		"ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"};

// Select global list by index (primary/secondary)
UniprotList * getGlobalLists(int passPrimary) {
	if (passPrimary) return uniprotPriEntryLists;
	else return uniprotSecEntryLists;
}

// Get global list count by index (primary/secondary)
int * getGlobalCount(int passPrimary) {
	if (passPrimary) return &uniprotPriListCount;
	else return &uniprotSecListCount;
}


void prepareUniprotReport(int passType, int passPrimary, UniprotList * passLists, CURLBuffer * passBuffer, const char * passProxy) {
	UniprotList * globalLists;

	// Aggregate and sort lists by value
	prepareUniprotLists(passLists, passPrimary);

	// Do not process if no results
	globalLists = getGlobalLists(passPrimary);
	if (globalLists[0].entryCount == 0) return;

	// Report specific preparation
	if (passType == OUTPUT_TYPE_UNIPROT_FULL) {
		// Submit entries to UniProt and retrieve full information
		retrieveUniprotOnline(passLists + UNIPROT_LIST_FULL, passBuffer, passProxy);
		joinOnlineLists(passLists + UNIPROT_LIST_FULL, passBuffer->buffer);
	}

	// Sort aggregated lists by count
	qsort(passLists[UNIPROT_LIST_FULL].entries, passLists[UNIPROT_LIST_FULL].entryCount, sizeof(UniprotEntry), uniprotEntryCompareID);
	qsort(passLists[UNIPROT_LIST_GENES].entries, passLists[UNIPROT_LIST_GENES].entryCount, sizeof(UniprotEntry), uniprotEntryCompareGene);
	qsort(passLists[UNIPROT_LIST_ORGANISM].entries, passLists[UNIPROT_LIST_ORGANISM].entryCount, sizeof(UniprotEntry), uniprotEntryCompareOrganism);
}

void renderUniprotReport(int passType, int passPrimary, FILE * passStream, const char * passProxy) {
	UniprotList * globalLists;
	UniprotList uniprotLists[3];
	CURLBuffer tempBuffer;
    char commonHeader[] = "Count\tAbundance\tQuality (Avg)\tQuality (Max)";

	// Prepare data
	prepareUniprotReport(passType, passPrimary, uniprotLists, &tempBuffer, passProxy);

	// Report no data
	globalLists = getGlobalLists(passPrimary);
	if (globalLists[0].entryCount == 0) {
		fprintf(passStream, "No entries to report\n");
		return;
	}

	// Render requested report
	switch (passType) {
	case OUTPUT_TYPE_UNIPROT_SIMPLE:
		fprintf(passStream, "%s\tUniProtKB\n", commonHeader);
		renderUniprotEntries(uniprotLists + UNIPROT_LIST_FULL, UNIPROT_LIST_FULL, passStream);
		fprintf(passStream, "\n\n%s\tGene\n", commonHeader);
		renderUniprotEntries(uniprotLists + UNIPROT_LIST_GENES, UNIPROT_LIST_GENES, passStream);
		fprintf(passStream, "\n\n%s\tOrganism\n", commonHeader);
		renderUniprotEntries(uniprotLists + UNIPROT_LIST_ORGANISM, UNIPROT_LIST_ORGANISM, passStream);

		break;

	case OUTPUT_TYPE_UNIPROT_FULL:
		fprintf(passStream, "%s\tUniProtKB\tID\tOrganism\tProtein Names\tGenes\tPathway\tFeatures\tGene Ontology\tReviewed\tExistence\tComments\tCross Reference (KEGG)\tCross Reference (GeneID)\tCross Reference (PATRIC)\tCross Reference(EnsemblBacteria)\tTaxonomic ID\tLineage\n", commonHeader);
		renderUniprotEntries(uniprotLists + UNIPROT_LIST_FULL, UNIPROT_LIST_FULL, passStream);
		freeCURLBuffer(&tempBuffer);

		break;
	}

	cleanUniprotLists(uniprotLists, passPrimary);
}

// Some UniProt references need to be reformatted for compatibility with reporting
int cleanUniprotReference(int passReference, const char * passBase) {
	char * tempName;
	int indexed;
	IndexHeader newHeader;
	FILE * proHandle;

	logMessage(__func__, LOG_LEVEL_MESSAGE, "Cleaning UniProt reference...\n");

	// Check if reference has already been indexed
	tempName = malloc(strlen(passBase) + 5);
	sprintf(tempName, "%s.ann", passBase);
	indexed = (access(tempName, F_OK) != -1);

	// Begin reference-specific preparation
	switch (passReference) {
	case UNIPROT_REFERENCE_UNIREF90:
		cleanUniprotReferenceUniref(passBase, 0);
		if (indexed) cleanUniprotReferenceUniref(tempName, 1);
	}

	// If indexed, also fix protein header
	if (indexed) {
		sprintf(tempName, "%s.pro", passBase);
		proHandle = err_xopen_core(__func__, tempName, "w");

		newHeader.multiFrame = 1;
		newHeader.nucleotide = 0;
		newHeader.referenceType = passReference;
		writeIndexHeader(proHandle, newHeader);

		fclose(proHandle);
	}

	free(tempName);

	return indexed;
}

// Rewrite headers with representative organism as UniProtKB ID for UniRef references.
void cleanUniprotReferenceUniref(const char * passName, int passANN) {
	gzFile srcHandle, dstHandle;
	FILE * dstANNHandle;
	int searchIdx, searchIdx2, lineIdx;
	char * newName;
	char srcBuffer[4096], dstBuffer[4096];
	char target[] = "RepID=";

	// Prepare and open reference and temporary destination
	newName = malloc(strlen(passName) + 5);
	sprintf(newName, "%s.tmp", passName);

	srcHandle = xzopen(passName, "r");
	if (!passANN) dstHandle = xzopen(newName, "w");
	else dstANNHandle = err_xopen_core(__func__, newName, "w");
	lineIdx = 0;

	while (gzgets(srcHandle, srcBuffer, 4096)) {
		// For FASTA files, must start with ">".  For ANN, must be an odd line
		if (((passANN == 0) && (srcBuffer[0] == '>')) ||
				((passANN == 1) && (lineIdx % 2 == 1))) {
			// Alter header - strip CRLF
			for (searchIdx = strlen(srcBuffer) - 1 ; searchIdx >= 0 ; searchIdx--) {
				if ((srcBuffer[searchIdx] == 0x0A) || (srcBuffer[searchIdx] == 0x0D)) {
					srcBuffer[searchIdx] = 0x00;
				}
				else break;
			}

			// Find RepID
			for (searchIdx = 0 ; searchIdx < strlen(srcBuffer) - strlen(target) ; searchIdx++) {
				for (searchIdx2 = 0 ; searchIdx2 < strlen(target) ; searchIdx2++ ) {
					if (srcBuffer[searchIdx + searchIdx2] != target[searchIdx2]) break;
				}

				if (searchIdx2 == strlen(target)) {
					// RepID found, write new header
					if (!passANN) sprintf(dstBuffer, ">%s %s\n", srcBuffer + searchIdx + strlen(target), srcBuffer + 1);
					else sprintf(dstBuffer, "0 %s %s\n", srcBuffer + searchIdx + strlen(target), srcBuffer + 2);
					break;
				}
			}

			if (!passANN) gzwrite(dstHandle, dstBuffer, strlen(dstBuffer));
			else fwrite(dstBuffer, 1, strlen(dstBuffer), dstANNHandle);
		}
		else {
			// Copy sequence
			if (!passANN) gzwrite(dstHandle, srcBuffer, strlen(srcBuffer));
			else fwrite(srcBuffer, 1, strlen(srcBuffer), dstANNHandle);
		}

		lineIdx++;
	}

	err_gzclose(srcHandle);
	if (!passANN) err_gzclose(dstHandle);
	else fclose(dstANNHandle);

	// Delete original file and replace with new one
	unlink(passName);
	rename(newName, passName);
}

// Download the requested UniProt reference (sprot/trembl/uniref90)
const char * downloadUniprotReference(int passReference, const char * passProxy) {
	CURL * curlHandle;
	CURLcode curlResult;
	FILE * fileHandle;
	const char * retFile;

	curlHandle = curl_easy_init();
	fileHandle = err_xopen_core(__func__, downloadNames[passReference], "w");
	retFile = downloadNames[passReference];

	curl_easy_setopt(curlHandle, CURLOPT_URL, downloadURLs[passReference]);
	curl_easy_setopt(curlHandle, CURLOPT_WRITEDATA, fileHandle);
    if (passProxy) curl_easy_setopt(curlHandle, CURLOPT_PROXY, passProxy);

	logMessage(__func__, LOG_LEVEL_MESSAGE, "Downloading %s...\n", downloadURLs[passReference]);
	curlResult = curl_easy_perform(curlHandle);

	if (curlResult != CURLE_OK) {
		logMessage(__func__, LOG_LEVEL_ERROR, "%s\n", curl_easy_strerror(curlResult));
		unlink(downloadNames[passReference]);
		retFile = "";
	}

	curl_easy_cleanup(curlHandle);
	fclose(fileHandle);

	return retFile;
}

void retrieveUniprotOnline(UniprotList * passList, CURLBuffer * retBuffer, const char * passProxy) {
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
	queryCount = (passList->entryCount < UNIPROT_MAX_SUBMIT) ? passList->entryCount : UNIPROT_MAX_SUBMIT;

	for (entryIdx = 0 ; entryIdx < passList->entryCount ; ) {
        // Preparation - build and sanitize query string outside of error loop 
        queryString[0] = 0;
        for (queryIdx = 0 ; (queryIdx < queryCount) && (entryIdx < passList->entryCount) ; entryIdx++) {
            for (parseIdx = 0 ; parseIdx < strlen(passList->entries[entryIdx].id) ; parseIdx++) {
                sprintf(queryString, "%s%s ", queryString, passList->entries[entryIdx].id);
                queryIdx++;
                break;
            }
        }

        httpString = curl_easy_escape(curlHandle, queryString, 0);
        sprintf(queryString, "uploadQuery=%s&format=job&from=ACC+ID&to=ACC&landingPage=false", httpString);
        curl_free(httpString);

		// Restart a limited number of times if errors encountered
		for (errorIdx = 0 ; errorIdx < UNIPROT_MAX_ERROR ; errorIdx++) {
            // Stage 1 - Submit query for processing
            logMessage(__func__, LOG_LEVEL_MESSAGE, "Submitting %d of %d entries to UniProt...\n", entryIdx, passList->entryCount);

			curl_easy_setopt(curlHandle, CURLOPT_URL, "https://www.uniprot.org/uploadlists/");
			curl_easy_setopt(curlHandle, CURLOPT_POSTFIELDS, queryString);
			curl_easy_setopt(curlHandle, CURLOPT_FOLLOWLOCATION, 1L);
			curl_easy_setopt(curlHandle, CURLOPT_WRITEFUNCTION, receiveUniprotOutput);
			curl_easy_setopt(curlHandle, CURLOPT_WRITEDATA, &tempBuffer);
            if (passProxy) curl_easy_setopt(curlHandle, CURLOPT_PROXY, passProxy);

			resetCURLBuffer(&tempBuffer);
			curlResult = curl_easy_perform(curlHandle);

			if (curlResult != CURLE_OK) {
				logMessage(__func__, LOG_LEVEL_ERROR, "CURL: %s (retrying)\n", curl_easy_strerror(curlResult));
				continue;
			}

			// Stage 2 - Wait for results
			if (tempBuffer.size > 50) {
				logMessage(__func__, LOG_LEVEL_ERROR, "Received unexpected job ID size (retrying)\n");
				continue;
			}
			sprintf(jobID, "%s", tempBuffer.buffer);
			sprintf(queryString, "https://www.uniprot.org/jobs/%s.stat", jobID);

			resetCURLBuffer(&tempBuffer);
			curl_easy_setopt(curlHandle, CURLOPT_URL, queryString);

			while (strcmp(tempBuffer.buffer, "COMPLETED") != 0) {
				sleep(1);
				resetCURLBuffer(&tempBuffer);
				curlResult = curl_easy_perform(curlHandle);
				if (curlResult != CURLE_OK) break;
                if (tempBuffer.size > 9 ) break;
			}

            if (tempBuffer.size > 9) {
                logMessage(__func__, LOG_LEVEL_ERROR, "Received unexpected job response sizei (retrying)\n");
                continue;
            }
			if (curlResult != CURLE_OK) {
				logMessage(__func__, LOG_LEVEL_ERROR, "CURL: %s (retrying)\n", curl_easy_strerror(curlResult));
				continue;
			}

			// Stage 3 - Retrieve Results
			sprintf(queryString, "query=job:%s&format=tab&columns=entry%%20name,id,organism,protein%%20names,genes,pathway,features,go,reviewed,existence,comments,database(KEGG),database(GeneID),database(PATRIC),database(EnsemblBacteria),organism-id,lineage(all)", jobID);
			curl_easy_setopt(curlHandle, CURLOPT_URL, "https://www.uniprot.org/uniprot/");
			curl_easy_setopt(curlHandle, CURLOPT_WRITEDATA, retBuffer);

			curlResult = curl_easy_perform(curlHandle);

			if (curlResult != CURLE_OK) {
				logMessage(__func__, LOG_LEVEL_ERROR, "CURL: %s\n", curl_easy_strerror(curlResult));
				continue;
			}

			break;
		}

        // Check if download failed max number of times
        if (errorIdx == UNIPROT_MAX_ERROR) {
            logMessage(__func__, LOG_LEVEL_ERROR, "Download failed after multiple retries.  Please check Internet connection and finalize with PALADIN-plugins\n");
            break;
        }
	}

	curl_easy_cleanup(curlHandle);
	freeCURLBuffer(&tempBuffer);
}

void renderUniprotEntries(UniprotList * passList, int passType, FILE * passStream) {
	int entryIdx, occurTotal;
	float occurPercent, avgQuality;
    char commonFields[] = "%d\t%.5f\t%.5f\t%d\t%s\n";

	// Count total occurrences for percentages
	for (entryIdx = 0, occurTotal = 0 ; entryIdx < passList->entryCount ; entryIdx++) {
		occurTotal += passList->entries[entryIdx].numOccurrence;
	}

	// Render fields
	for (entryIdx = 0 ; entryIdx < passList->entryCount ; entryIdx++) {
		occurPercent = (float) passList->entries[entryIdx].numOccurrence / (float) occurTotal * 100;
		avgQuality = (float) passList->entries[entryIdx].totalQuality / (float) passList->entries[entryIdx].numOccurrence;

		switch(passType) {
		    case UNIPROT_LIST_FULL:
		    	fprintf(passStream, commonFields, passList->entries[entryIdx].numOccurrence, occurPercent, avgQuality, passList->entries[entryIdx].maxQuality,  passList->entries[entryIdx].id);
		    	break;
		    case UNIPROT_LIST_GENES:
			    fprintf(passStream, commonFields, passList->entries[entryIdx].numOccurrence, occurPercent, avgQuality, passList->entries[entryIdx].maxQuality, passList->entries[entryIdx].gene);
			    break;
		    case UNIPROT_LIST_ORGANISM:
			    fprintf(passStream, commonFields, passList->entries[entryIdx].numOccurrence, occurPercent, avgQuality, passList->entries[entryIdx].maxQuality, passList->entries[entryIdx].organism);
			    break;
		}
	}
}

void renderNumberAligned(const mem_opt_t * passOptions) {
	int listIdx, successTotal, alignTotal;

	// Primary alignments
	for (successTotal = 0, alignTotal = 0, listIdx = 0 ; listIdx < uniprotPriListCount ; listIdx++) {
		successTotal += uniprotPriEntryLists[listIdx].entryCount;
		alignTotal += uniprotPriEntryLists[listIdx].entryCount + uniprotPriEntryLists[listIdx].unalignedCount;
	}

	// Secondary alignments (if requested)
	if (passOptions->flag & MEM_F_ALL) {
		for (listIdx = 0 ; listIdx < uniprotSecListCount ; listIdx++) {
			successTotal += uniprotSecEntryLists[listIdx].entryCount;
			alignTotal += uniprotSecEntryLists[listIdx].entryCount;
		}
	}

	if (alignTotal == 0) {
		logMessage(__func__, LOG_LEVEL_MESSAGE, "No detected ORF sequences, no alignment performed\n");
	}
	else {
		logMessage(__func__, LOG_LEVEL_MESSAGE, "Aligned %d out of %d total detected ORF sequences (%.2f%%)\n", successTotal, alignTotal, (float) successTotal / (float) alignTotal * 100);
	}
}

int addUniprotList(worker_t * passWorker, int passSize, int passFull) {
	int entryIdx, alnIdx, addPriIdx, addSecIdx, parseIdx;
	int refID, alignType, primaryCount, totalAlign, entryQuality;
	UniprotList * globalLists;
	int * globalCount, * currentIdx;
	char * uniprotEntry;

	// Create lists
	uniprotPriEntryLists = realloc(uniprotPriEntryLists, (uniprotPriListCount + 1) * sizeof(UniprotList));
	memset(uniprotPriEntryLists + uniprotPriListCount, 0, sizeof(UniprotList));
	uniprotSecEntryLists = realloc(uniprotSecEntryLists, (uniprotSecListCount + 1) * sizeof(UniprotList));
	memset(uniprotSecEntryLists + uniprotSecListCount, 0, sizeof(UniprotList));

	// Calculate potential total list size
	for (entryIdx = 0, totalAlign = 0 ; entryIdx < passSize ; entryIdx++) {
		// Only add active sequences
		if (!passWorker->regs[entryIdx].active) continue;

		for (alnIdx = 0, primaryCount = 0 ; alnIdx < passWorker->regs[entryIdx].n ; alnIdx++) {
			switch (alignType = getAlignmentType(passWorker, entryIdx, alnIdx)) {
			case MEM_ALIGN_PRIMARY:
				uniprotPriEntryLists[uniprotPriListCount].entryCount++;
				// Count non-linear as total alignments
				if (primaryCount++) totalAlign++;
				break;
			case MEM_ALIGN_SECONDARY:
				uniprotSecEntryLists[uniprotSecListCount].entryCount++; break;
			}
		}	

		totalAlign++;
	}

	// Set counts and allocate entries
	uniprotPriEntryLists[uniprotPriListCount].unalignedCount = totalAlign - uniprotPriEntryLists[uniprotPriListCount].entryCount;
	uniprotSecEntryLists[uniprotSecListCount].unalignedCount = totalAlign - uniprotSecEntryLists[uniprotSecListCount].entryCount;;
	uniprotPriEntryLists[uniprotPriListCount].entries = calloc(uniprotPriEntryLists[uniprotPriListCount].entryCount, sizeof(UniprotEntry));
	uniprotSecEntryLists[uniprotSecListCount].entries = calloc(uniprotSecEntryLists[uniprotSecListCount].entryCount, sizeof(UniprotEntry));

	// Populate list
	for (addPriIdx = 0, addSecIdx = 0, entryIdx = 0 ; entryIdx < passSize && passFull ; entryIdx++) {
		// Only add active sequences
		if (!passWorker->regs[entryIdx].active) continue;

		for (alnIdx = 0 ; alnIdx < passWorker->regs[entryIdx].n ; alnIdx++) {
			// Only add successful alignments
			if ((alignType = getAlignmentType(passWorker, entryIdx, alnIdx)) < MEM_ALIGN_PRIMARY) continue;

			globalLists = getGlobalLists(alignType == MEM_ALIGN_PRIMARY);
			globalCount = getGlobalCount(alignType == MEM_ALIGN_PRIMARY);
			currentIdx = (alignType == MEM_ALIGN_PRIMARY) ? &addPriIdx : &addSecIdx;

			// Extract ID and entry
			refID = passWorker->regs[entryIdx].a[alnIdx].rid;
			uniprotEntry = passWorker->bns->anns[refID].name;

			// Strip sequence and frame info for nucleotide references
			if (passWorker->opt->indexInfo.nucleotide) {
				for (parseIdx = 2 ; (parseIdx > 0) && (*uniprotEntry != 0) ; uniprotEntry++) {
					if (*uniprotEntry == ':') parseIdx--;
				}
			}

			// Strip initial IDs (if present) and description
			for (parseIdx = 0 ; (uniprotEntry[parseIdx] != 0) && (uniprotEntry[parseIdx] != ' ') ; parseIdx++) {
				if (uniprotEntry[parseIdx] == '|') {
					uniprotEntry += parseIdx + 1;
					parseIdx = 0;
				}
			}

			uniprotEntry[parseIdx] = 0;	

			// Full ID and quality
			globalLists[*globalCount].entries[*currentIdx].id = malloc(strlen(uniprotEntry) + 1);
			sprintf(globalLists[*globalCount].entries[*currentIdx].id, "%s", uniprotEntry);
			globalLists[*globalCount].entries[*currentIdx].numOccurrence = 1;

            entryQuality = passWorker->regs[entryIdx].a[alnIdx].mapq;
			globalLists[*globalCount].entries[*currentIdx].totalQuality = entryQuality;
            if (entryQuality > globalLists[*globalCount].entries[*currentIdx].maxQuality) {
                globalLists[*globalCount].entries[*currentIdx].maxQuality = entryQuality;
            }

			// Gene/organism
			for (parseIdx = 0 ; parseIdx < strlen(uniprotEntry) ; parseIdx++) {
				if (*(uniprotEntry + parseIdx) == '_') {
					globalLists[*globalCount].entries[*currentIdx].gene = malloc(parseIdx + 1);
					sprintf(globalLists[*globalCount].entries[*currentIdx].gene, "%.*s", parseIdx, uniprotEntry);
					globalLists[*globalCount].entries[*currentIdx].organism = malloc(strlen(uniprotEntry + parseIdx) + 1);
					sprintf(globalLists[*globalCount].entries[*currentIdx].organism, "%s", uniprotEntry + parseIdx + 1);
					parseIdx = -1;
					break;
				}
			}

			// If underscore missing, we may be dealing with clustered ID with deleted representative
			if (parseIdx > -1) {
				globalLists[*globalCount].entries[*currentIdx].gene = malloc(strlen(uniprotEntry) + 1);
				sprintf(globalLists[*globalCount].entries[*currentIdx].gene, "%s", uniprotEntry);
				globalLists[*globalCount].entries[*currentIdx].organism = malloc(8);
				sprintf(globalLists[*globalCount].entries[*currentIdx].organism, "Unknown");
			}

			(*currentIdx)++;
		}
	}

	uniprotPriListCount++;
	uniprotSecListCount++;

	return uniprotPriListCount - 1;
}

void cleanUniprotLists(UniprotList * passLists, int passPrimary) {
	UniprotList * globalLists;
	int listIdx, entryIdx;

	globalLists = getGlobalLists(passPrimary);

	// Global list 0 contains pointers to all allocated strings, so delete from there. Do not repeat with other lists
	for (entryIdx = 0 ; entryIdx < globalLists[0].entryCount ; entryIdx++) {
		free(globalLists[0].entries[entryIdx].id);
		free(globalLists[0].entries[entryIdx].gene);
		free(globalLists[0].entries[entryIdx].organism);
	}

	free(globalLists[0].entries);
	free(globalLists);

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


void prepareUniprotLists(UniprotList * retLists, int passPrimary) {
	UniprotList * globalLists;
	int listIdx, entryIdx, localIdx;
	int maxEntries, totalSize;

	globalLists = getGlobalLists(passPrimary);

	// Calculate maximum number of UniProt entries
	for (maxEntries = 0, listIdx = 0 ; listIdx < *getGlobalCount(passPrimary) ; listIdx++) {
		maxEntries += globalLists[listIdx].entryCount;
	}

	logMessage(__func__, LOG_LEVEL_MESSAGE, "Aggregating %d entries for UniProt report\n", maxEntries);

	// Stop processing if no entries, but ensure a list exists for easier post-processing
	if (maxEntries == 0) {
		if (globalLists == NULL) {
			if (passPrimary) uniprotPriEntryLists = malloc(sizeof(UniprotList));
			else uniprotSecEntryLists = malloc(sizeof(UniprotList));

			globalLists = getGlobalLists(passPrimary);
			globalLists[0].entries = NULL;
			globalLists[0].entryCount = 0;
			globalLists[0].unalignedCount = 0;
			*getGlobalCount(passPrimary) = 1;
		}

		return;
	}

	// Join each pipeline's individual list into one master list (at first list)
	for (listIdx = 0, totalSize = 0 ; listIdx < *getGlobalCount(passPrimary) ; listIdx++) {
		totalSize += globalLists[listIdx].entryCount;
	}

	globalLists[0].entries = realloc(globalLists[0].entries, sizeof(UniprotEntry) * totalSize);

	for (listIdx = 1, entryIdx = globalLists[0].entryCount ; listIdx < *getGlobalCount(passPrimary) ; listIdx++) {
		memcpy(globalLists[0].entries + entryIdx, globalLists[listIdx].entries, sizeof(UniprotEntry) * globalLists[listIdx].entryCount);
		entryIdx += globalLists[listIdx].entryCount;
		free(globalLists[listIdx].entries);
		globalLists[listIdx].entryCount = 0;
	}

	globalLists[0].entryCount = totalSize;

	// Aggregate each local list
	for (localIdx = 0 ; localIdx < 3 ; localIdx++) {
		aggregateUniprotList(retLists + localIdx, localIdx, passPrimary);
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
	for (entryIdx = 0, lineIdx = 1 ; (entryIdx < retList->entryCount) && (lineIdx < lineCount)  ; ) {
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

void aggregateUniprotList(UniprotList * retList, int passListType, int passPrimary) {
	UniprotList * globalLists;
	int entryIdx;
	int memberOffset;
	char * srcBase, * dstBase;

	globalLists = getGlobalLists(passPrimary);

	// First sort full list
	switch (passListType) {
	case UNIPROT_LIST_FULL:
		qsort(globalLists[0].entries, globalLists[0].entryCount, sizeof(UniprotEntry), uniprotEntryCompareID);
		memberOffset = offsetof(UniprotEntry, id);
		break;
	case UNIPROT_LIST_GENES:
		qsort(globalLists[0].entries, globalLists[0].entryCount, sizeof(UniprotEntry), uniprotEntryCompareGene);
		memberOffset = offsetof(UniprotEntry, gene);
		break;
	case UNIPROT_LIST_ORGANISM:
		qsort(globalLists[0].entries, globalLists[0].entryCount, sizeof(UniprotEntry), uniprotEntryCompareOrganism);
		memberOffset = offsetof(UniprotEntry, organism);
		break;
	default: memberOffset = 0; break;
	}

	// Count entry occurrences and aggregate into return list
	retList->entries = calloc(globalLists[0].entryCount, sizeof(UniprotEntry));
	retList->entries[0].id = globalLists[0].entries[0].id;
	retList->entries[0].gene = globalLists[0].entries[0].gene;
	retList->entries[0].organism = globalLists[0].entries[0].organism;
	retList->entryCount = 0;

	for (entryIdx = 0 ; entryIdx < globalLists[0].entryCount ; entryIdx++) {
		srcBase = (char *) (globalLists[0].entries + entryIdx);
		dstBase = (char *) (retList->entries + retList->entryCount);

		// If current global entry doesn't match local entry, create new local entry and copy pointer to relevant string info
		if (strcmp(*((char * *) (srcBase + memberOffset)), *((char * *) (dstBase + memberOffset))) != 0) {
			dstBase = (char *) (retList->entries + ++retList->entryCount);
			*((char * *) (dstBase + memberOffset)) = *((char * *) (srcBase + memberOffset));
		}

		retList->entries[retList->entryCount].numOccurrence++;
		retList->entries[retList->entryCount].totalQuality += (globalLists[0].entries + entryIdx)->totalQuality;
        if ((globalLists[0].entries + entryIdx)->maxQuality > retList->entries[retList->entryCount].maxQuality) {
            retList->entries[retList->entryCount].maxQuality = (globalLists[0].entries + entryIdx)->maxQuality;
        }
	}

	retList->entryCount++;
}


int uniprotEntryCompareID (const void * passEntry1, const void * passEntry2) {
	int compVal = 0;

	if (((compVal = ((UniprotEntry *)passEntry2)->numOccurrence - ((UniprotEntry *)passEntry1)->numOccurrence)) != 0) {
		return compVal;
	}

	if ((compVal = strcmp(((UniprotEntry *)passEntry1)->id, ((UniprotEntry *)passEntry2)->id)) != 0) {
		return compVal;
	}

	return compVal;
}

int uniprotEntryCompareGene (const void * passEntry1, const void * passEntry2) {
	int compVal = 0;

	if (((compVal = ((UniprotEntry *)passEntry2)->numOccurrence - ((UniprotEntry *)passEntry1)->numOccurrence)) != 0) {
		return compVal;
	}

	if ((compVal = strcmp(((UniprotEntry *)passEntry1)->gene, ((UniprotEntry *)passEntry2)->gene)) != 0) {
		return compVal;
	}

	return compVal;
}

int uniprotEntryCompareOrganism (const void * passEntry1, const void * passEntry2) {
	int compVal = 0;

	if (((compVal = ((UniprotEntry *)passEntry2)->numOccurrence - ((UniprotEntry *)passEntry1)->numOccurrence)) != 0) {
		return compVal;
	}

	if ((compVal = strcmp(((UniprotEntry *)passEntry1)->organism, ((UniprotEntry *)passEntry2)->organism)) != 0) {
		return compVal;
	}

	return compVal;
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
