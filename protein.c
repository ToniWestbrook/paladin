/* Contact: Toni Westbrook <anthonyw@wildcats.unh.edu> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include "protein.h"
#include "bntseq.h"
#include "utils.h"
#include "bwt.h"
#include "main.h"

#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif
extern unsigned char nst_nt4_table[256];

// Codon encoded as 6 bit value, MSB as left-most nucleotide in 3-mer
unsigned char codon_aa_hash[64] = {
	'K', 'N', 'K', 'N', // AA?
	'T', 'T', 'T', 'T', // AC?
	'R', 'S', 'R', 'S', // AG?
	'I', 'I', 'M', 'I', // AT?

	'Q', 'H', 'Q', 'H', // CA?
	'P', 'P', 'P', 'P', // CC?
	'R', 'R', 'R', 'R', // CG?
	'L', 'L', 'L', 'L', // CT?

	'E', 'D', 'E', 'D', // GA?
	'A', 'A', 'A', 'A', // GC?
	'G', 'G', 'G', 'G', // GG?
	'V', 'V', 'V', 'V', // GT?

	'*', 'Y', '*', 'Y', // TA?
	'S', 'S', 'S', 'S', // TC?
	'*', 'C', 'W', 'C', // TG?
	'L', 'F', 'L', 'F', // TT?
};

// Debug/Test function to lookup corresponding CDS to start read location
void testReadHeader(char * passOrig, char * retNew) {
	FILE * gffHandle;
	int parseIdx, colonCount;
	unsigned long searchIdx, line;
	struct CDS currentCDS;

	retNew[0] = 0;

	// Find start of number
	for (parseIdx = 0, colonCount = 0 ; parseIdx < strlen(passOrig) ; parseIdx++) {
		if (passOrig[parseIdx] == ':') colonCount++;
		if (colonCount == 10) break;
	}

	// Found number
	parseIdx++;
	if (parseIdx < strlen(passOrig)) {
		searchIdx = atol(passOrig + parseIdx);

		gffHandle = fopen("/thomas1/mcbs913/anthony/mapping/degen/panda/EscherichiaColiStrK-12SubstrMG1655/ecoli.gff", "r");
		while (getNextCDS(gffHandle, &currentCDS, &line)) {
			if ((searchIdx >= currentCDS.startIdx) && (searchIdx < currentCDS.endIdx)) {
				sprintf(retNew, "%s", currentCDS.description);
				break;
			}
		}
		fclose(gffHandle);
	}
}

// Write header for index pro file
void writeIndexHeader(FILE * passFilePtr, int passProtein, int passMulti) {
	fprintf(passFilePtr, ">NT=%d:MF=%d:VER=%s\n", passProtein, passMulti, PACKAGE_VERSION);
}

// Get header info from index pro file
char getIndexHeader(char * passFile) {
	FILE * filePtr;
	int readNT, readMF;
	char retValue;
	
	// Open file and read header information
	filePtr = fopen(passFile, "r");
	
    if (filePtr == NULL) {
    	fprintf(stderr, "[%s] fail to open file '%s' : %s\n", __func__, passFile, errno ? strerror(errno) : "Out of memory");
    	exit(EXIT_FAILURE);
    }


	if (fscanf(filePtr, ">NT=%d:MF=%d", &readNT, &readMF) != 2) {
    	fprintf(stderr, "[%s] fail to parse file '%s'\n", __func__, passFile);
    	exit(EXIT_FAILURE);
	}
	fclose(filePtr);
	
	retValue = 0;
	if (readNT) retValue |= INDEX_FLAG_NT;
	if (readMF) retValue |= INDEX_FLAG_MF;

	return retValue;
}


// Construct 6-bit codon from 3 ASCII characters
unsigned char encodeCodon(char * passSequence, int passStrand) {
	unsigned char retCodon;
	int relIdx, absIdx;

	retCodon = 0;

	for (relIdx = 0 ; relIdx < 3 ; relIdx++) {
		// Calculate absolute index from relative
		absIdx = relIdx;
		if (passStrand < 0) absIdx = absIdx * -1;
 
		// If this codon has any ambiguous nucleotides, signal with 0xFF
		if (nst_nt4_table[(int) passSequence[absIdx]] > 3) {
			return 0xFF;
		}

		retCodon = (retCodon << 2) | nst_nt4_table[(int) passSequence[absIdx]];
	}

	// Complement if reverse strand
	if (passStrand < 0) retCodon ^= 0x3F;

	return retCodon;
}

// Give a nucleotide sequence and a CDS entry (containing alignment info), return the corresponding protein sequence
int convertToAA(char * passSequence, struct CDS * passCDS, char ** retSequence, unsigned long * retSize) {
	unsigned long nucIdx, aaIdx;
	unsigned long seqLen, seqStart;
	unsigned char currentCodon;

	// Calculate sequence length and loop control from CDS entry
	seqLen = passCDS->endIdx + 1 - passCDS->startIdx - passCDS->phase;
	seqStart = (passCDS->strand == 1) ? passCDS->startIdx : passCDS->endIdx;

	// Reserve memory for amino acid sequence
	*retSize = (seqLen + 2) / 3;
	*retSequence = malloc(*retSize);

	// Iterate through all codons in the nucleotide sequence
	for (nucIdx = 0, aaIdx = 0 ; nucIdx < seqLen ; nucIdx += 3, aaIdx++) {
		// If we've exhausted nucleotides premature to a full frame, assign ambiguous AA and exit
		if (nucIdx + 3 > seqLen) {
			(*retSequence)[aaIdx] = 'X';
			break;
		}

		// Encode codon
		currentCodon = encodeCodon(passSequence + seqStart + (passCDS->strand * nucIdx), passCDS->strand);
		
		// If ambiguous nucleotides are present, assign ambiguous AA
		if (currentCodon == 0xFF) {
			(*retSequence)[aaIdx] = 'X';
		}

		// Hash to amino acid IUPAC
		(*retSequence)[aaIdx] = codon_aa_hash[currentCodon];
	}
	
	return 0;
}

// Calculate last full position of ORF 
long getLastORFPos(long passLength, int passFrame) {
	long retPos;

	if (passFrame < 3) {
		// Forward frame
		retPos = passLength - (passLength % 3) - 1 + (passFrame % 3);
		if (retPos >= passLength) retPos -= 3;
	}
	else {
		// Reverse frame
		retPos = (passLength % 3) - (passFrame % 3);
		if (retPos < 0) retPos += 3;
	}

	return retPos;
}


// Add new entry to dynamic history array
void addORFHistory(long * passHistory[2][6], long passHistorySize[6], unsigned long passIdx ) {
	int copyIdx;
	long * tempHistory;

	if (passHistorySize[passIdx] == 0) {
		// New array
		passHistory[0][passIdx] = malloc(sizeof(long));
		passHistory[1][passIdx] = malloc(sizeof(long));
	}
	else {
		// If array has existing items, resize, copy, and free old
		for (copyIdx = 0 ; copyIdx < 2 ; copyIdx++) {
			tempHistory = passHistory[copyIdx][passIdx];
			passHistory[copyIdx][passIdx] = malloc((passHistorySize[passIdx] + 1) * sizeof(long));
			memcpy(passHistory[copyIdx][passIdx], tempHistory, passHistorySize[passIdx] * sizeof(long));
			free(tempHistory);
		}
	}

	// Update size
	passHistorySize[passIdx]++;
}

// Filter history array for ORF selection parameters, return corresponding CDS array
void compileORFHistory(long * passHistory[2][6], long passHistorySize[6], struct CDS * * retCDS, unsigned long * retCount) {
	unsigned long srcFrameIdx, srcHistoryIdx;
	unsigned long srcStart, srcEnd, totalSize;
	int relStart, validEntry;

	// Allocate potential CDS records per entry in history array
	totalSize = 0;
	for (srcFrameIdx = 0 ; srcFrameIdx < 6 ; srcFrameIdx++) totalSize += passHistorySize[srcFrameIdx];

	*retCDS = calloc(totalSize, sizeof(struct CDS));
	*retCount = 0;
	relStart = -1;

	// Iterate through each entry in each reading frame
	for (srcFrameIdx = 0 ; srcFrameIdx < 6 ; srcFrameIdx++) {
		for (srcHistoryIdx = 0 ; srcHistoryIdx < passHistorySize[srcFrameIdx] ; srcHistoryIdx++) {
			// Swap start/stop for reverse direction
			if (passHistory[0][srcFrameIdx][srcHistoryIdx] < passHistory[1][srcFrameIdx][srcHistoryIdx]) {
				srcStart = passHistory[0][srcFrameIdx][srcHistoryIdx];
				srcEnd = passHistory[1][srcFrameIdx][srcHistoryIdx];
			}
			else {
				srcStart = passHistory[1][srcFrameIdx][srcHistoryIdx];
				srcEnd = passHistory[0][srcFrameIdx][srcHistoryIdx];
			}

			validEntry = 1;

			// Filter minimum ORF size
			//if (srcLen < passMinORF) continue;
			
			// Add CDS record if valid (current strategy marks all as valid)
			if (validEntry) {
				if (relStart < 0) relStart = srcFrameIdx;

				(*retCDS)[*retCount].startIdx = srcStart;
				(*retCDS)[*retCount].endIdx = srcEnd;
				(*retCDS)[*retCount].strand = (srcFrameIdx < 3) ? 1 : -1;
				(*retCDS)[*retCount].relFrame = srcFrameIdx - relStart;

				(*retCount)++;
			}
		}
	}
}

// Scan nucleotide sequence for all recognized ORFs, return as CDS array
int getSequenceORF(char * passSequence, unsigned long passLength, mem_opt_t * passOptions, struct CDS * * retCDS, unsigned long * retCount) {
	int frameIdx, strandDir, relFrame, searchIdx, highIdx;
	long * history[2][6], historySize[6], stopCounts[6], stopOrder[6];
	long gcCount, seqIdx, absIdx;
	unsigned char currentCodon;

	// Initialize values
	gcCount = 0;
	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		history[0][frameIdx] = NULL;
		history[1][frameIdx] = NULL;
		historySize[frameIdx] = 0;
		stopCounts[frameIdx] = 0;
		stopOrder[frameIdx] = 0;
	}

	// Collect statistics (stop codon counts and GC content)
	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		strandDir = (frameIdx < 3 ? 1 : -1);
		relFrame = frameIdx % 3;

		// Iterate through all codons in this read frame
		for (seqIdx = (frameIdx % 3) ; seqIdx + 3 <= passLength ; seqIdx += 3) {
			// Calculate absolute index from relative
			absIdx = seqIdx;
			if (strandDir == -1) absIdx = absIdx * -1 + passLength - 1;

			// Encode codon
			currentCodon = encodeCodon(passSequence + absIdx, strandDir);

			// Record occurrence of stop codon (TAA - 0x30, TAG - 0x32, TGA - 0x38)
			if ((currentCodon == 0x30) || (currentCodon == 0x32) || (currentCodon == 0x38)) {
				stopCounts[frameIdx]++;
				stopOrder[frameIdx]++;
			}

			// Detect GC content on frame 0
			if (frameIdx == 0) {
				if (((currentCodon & 0x03) == 0x01) || ((currentCodon & 0x03) == 0x02)) gcCount++;
				if (((currentCodon & 0x0C) == 0x04) || ((currentCodon & 0x0C) == 0x08)) gcCount++;
				if (((currentCodon & 0x30) == 0x10) || ((currentCodon & 0x30) == 0x20)) gcCount++;
			}
		}
	}

	// Use occurrence counts to generate frame order 
	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		highIdx = 0;

		for (searchIdx = 0 ; searchIdx < 6 ; searchIdx++) {
			if (stopOrder[searchIdx] > stopOrder[highIdx]) {
				highIdx = searchIdx;
			}
		}

		// Record order as negatives for inline efficiency, then convert after
		stopOrder[highIdx] = -5 + frameIdx;
	}

	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		stopOrder[frameIdx] *= -1;
	}

	// All modes iterate through all frames
	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		strandDir = (frameIdx < 3 ? 1 : -1);
		relFrame = frameIdx % 3;

		// Index is multiframe, stop at first ORF
		if (passOptions->indexFlag & INDEX_FLAG_MF) {
			if (stopCounts[frameIdx] == 0) {
				addORFHistory(history, historySize, frameIdx);
				history[0][frameIdx][historySize[frameIdx] - 1] = (strandDir == 1) ? relFrame : passLength - 1 - relFrame;
				history[1][frameIdx][historySize[frameIdx] - 1] = getLastORFPos(passLength, frameIdx);
				break;
			}
		}

		// Index is single frame, check for algorithm variant
		else {
			if (passOptions->proteinFlag & ALIGN_FLAG_BRUTEORF) {
				// Brute force ORF detection - encode all frames if ORF detected
				if (stopCounts[frameIdx] == 0) {
					for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
						strandDir = (frameIdx < 3 ? 1 : -1);
						relFrame = frameIdx % 3;

						addORFHistory(history, historySize, frameIdx);
						history[0][frameIdx][historySize[frameIdx] - 1] = (strandDir == 1) ? relFrame : passLength - 1 - relFrame;
						history[1][frameIdx][historySize[frameIdx] - 1] = getLastORFPos(passLength, frameIdx);
					}

					break;
				}
			}
			else {
				// Smart ORF detection - to be filled in
				if (stopCounts[frameIdx] <= 0) {
					if ((stopOrder[(frameIdx + 4) % 6] != 0) && (stopOrder[(frameIdx + 5) % 6] != 0)) {
						addORFHistory(history, historySize, frameIdx);
						history[0][frameIdx][historySize[frameIdx] - 1] = (strandDir == 1) ? relFrame : passLength - 1 - relFrame;
						history[1][frameIdx][historySize[frameIdx] - 1] = getLastORFPos(passLength, frameIdx);
					}
				}
			}
		}
	}

	// Filter, remove overlaps, create CDS array
	compileORFHistory(history, historySize, retCDS, retCount);

	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		free(history[0][frameIdx]);
		free(history[1][frameIdx]);
	}

	return 0;


		// 63%
		//if ((stopCounts[(frameIdx + 0) % 6] <= 1) &&
		//	((stopOrder[(frameIdx + 4) % 6] != 0) && (stopOrder[(frameIdx + 5) % 6] != 0))){


		// 50%
		//if (((stopOrder[(frameIdx + 0) % 6] == 5) && (stopOrder[(frameIdx + 3) % 6] == 4)) ||
		//	((stopOrder[(frameIdx + 0) % 6] == 4) && (stopOrder[(frameIdx + 3) % 6] == 5))){


		// First full ORF
//		if ((stopCounts[frameIdx] == 0)) {
/*			(stopOrder[(frameIdx + 0) % 6] <= 1) && (stopOrder[(frameIdx + 3) % 6] <= 1) &&
			(stopOrder[(frameIdx + 4) % 6] >= 4) && (stopOrder[(frameIdx + 5) % 6] >= 4) &&
			((stopOrder[(frameIdx + 1) % 6] == 2) || (stopOrder[(frameIdx + 1) % 6] == 3)) &&
			((stopOrder[(frameIdx + 2) % 6] == 2) || (stopOrder[(frameIdx + 2) % 6] == 3))) {*/
		//if ((stopCounts[(frameIdx + 0) % 6] <= 1) && (stopCounts[(frameIdx + 3) % 6] <= 1)) {

		//if ((stopCounts[frameIdx] < 1) &&
			//((stopOrder[frameIdx] == 0) || (stopOrder[frameIdx] == 1))){ //&&
			//((stopCounts[(frameIdx + 3) % 6] == 0) || (stopCounts[(frameIdx + 3) % 6] == 1))) {// &&
			//((stopCounts[(frameIdx + 2) % 6] == 2) || (stopCounts[(frameIdx + 2) % 6] == 3)) &&
			//(stopCounts[(frameIdx + 2) % 6] == 2) && (stopCounts[(frameIdx + 4) % 6] == 5)) {
			//((stopOrder[(frameIdx + 4) % 6] == 4) || (stopOrder[(frameIdx + 4) % 6] == 5))) {

	// Do something with this
	/*
	if (!stopFound) {
		addORFHistory(history, historySize, historyIdx);
		history[0][historyIdx][historySize[historyIdx] - 1] = (strandIdx == 0) ? frameIdx : passLength - 1 - frameIdx;
		history[1][historyIdx][historySize[historyIdx] - 1] = (strandIdx == 0) ? passLength - 4 + frameIdx : frameIdx;
		compileORFHistory(history, historySize, passMinORF, retCDS, retCount);
		return 0;
	}
	*/

}

// Iterate through GFF annotation file, return next available CDS entry
int getNextCDS(FILE * passFile, struct CDS * retCDS, unsigned long * retLine) {
	int fieldIdx, scanIdx;
	int readCount;
	char readBuffer[GFF_MAX_FIELD], field[9][GFF_MAX_FIELD];
	
	// Iterate through each annotation line until CDS or EOF
	while (fgets(readBuffer, GFF_MAX_FIELD, passFile)) {
		// Skip comments and blank lines
		(*retLine)++;
		if (readBuffer[0] == '#') continue;
		if ((readBuffer[0] == '\n') || (readBuffer[0] == '\r')) continue;

		// Parse all fields from GFF entry
		for (fieldIdx = 0, scanIdx = 0 ; fieldIdx < 9 ; fieldIdx++) {
			sscanf(readBuffer + scanIdx, "%[^\t]%n", field[fieldIdx], &readCount);
			scanIdx += readCount + 1;
		}

		// Check for CDS entry, assign to record and return
		if (strcmp(field[2], "CDS") == 0) {
			retCDS->startIdx = atol(field[3]) - 1;
			retCDS->endIdx = atol(field[4]) - 1;
			retCDS->strand = (field[6][0] == '+') ? 1 : -1;
			retCDS->phase = atoi(field[7]);

			// Fixup description
			sprintf(retCDS->description, "%s", field[8]);
			for (scanIdx = 0 ; scanIdx < strlen(retCDS->description) ; scanIdx++) {
				if (retCDS->description[scanIdx] == ' ') retCDS->description[scanIdx] = '_';
				if (retCDS->description[scanIdx] == 0x0A || retCDS->description[scanIdx] == 0x0D) retCDS->description[scanIdx] = 0;
			}

			return 1;
		}
	}

	return 0;
}

// Convert the given reference nucleotide FASTA and GFF file to a protein FASTA file (for all one/all frames)
int writeIndexProtein(const char * passPrefix, const char * passProName, const char * passAnnName, int passMulti) {
        struct CDS currentCDS;
        gzFile inputSeqPtr;
        FILE * inputAnnPtr, * outputPtr;
        char * outputBuffer;
        kseq_t * seq;
        unsigned long outputSize, currentLine, frameIdx;

        // Prepare file handles
        inputSeqPtr = xzopen(passPrefix, "r");
        inputAnnPtr = fopen(passAnnName, "r");
        outputPtr = fopen(passProName, "w");

        if (inputAnnPtr == NULL) {
        	fprintf(stderr, "[%s] fail to open file '%s' : %s\n", __func__, passAnnName, errno ? strerror(errno) : "Out of memory");
        	exit(EXIT_FAILURE);
        }

    	// Write index type header
    	writeIndexHeader(outputPtr, 1, passMulti);

        // Read in 1st sequence data
        seq = kseq_init(inputSeqPtr);
        kseq_read(seq);

        // Iterate through each CDS sequence in the annotation file
        currentLine = 0;
        while (getNextCDS(inputAnnPtr, &currentCDS, &currentLine)) {
        	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
				convertToAA(seq->seq.s, &currentCDS, &outputBuffer, &outputSize);

				// Line ID of CDS : Frame : Sequence Header
				fprintf(outputPtr, ">%lu:%lu:%s\n%.*s\n", currentLine, frameIdx, currentCDS.description, (int) outputSize, outputBuffer);
				free(outputBuffer);

                currentCDS.startIdx++;
				currentCDS.endIdx++;
				if (frameIdx == 2) currentCDS.strand *= -1;

				if (!passMulti) break;
			}
        }

        // Close files
        fflush(outputPtr);
        err_gzclose(inputSeqPtr);
        err_fclose(inputAnnPtr);
        err_fclose(outputPtr);
	kseq_destroy(seq);
        return 0;
}

int writeIndexCodingProtein(const char * passPrefix, const char * passProName, int passMulti) {
    struct CDS currentCDS;
    gzFile  inputSeqPtr;
    FILE * outputPtr;
    char * outputBuffer;
    kseq_t * seq;
    unsigned long outputSize, currentLine, frameIdx;

    // Prepare file handles
    inputSeqPtr = xzopen(passPrefix, "r");
    outputPtr = fopen(passProName, "w");

	// Write index type header
	writeIndexHeader(outputPtr, 1, passMulti);

    // Iterate through each sequence
    seq = kseq_init(inputSeqPtr);
    currentLine = 0;

    while (kseq_read(seq) > 0) {
            for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
            	if (frameIdx < 3) {
                    currentCDS.startIdx = (frameIdx % 3);
                    currentCDS.endIdx = getLastORFPos(seq->seq.l, frameIdx);
                    currentCDS.strand = 1;
                    currentCDS.phase = 0;
            	}
            	else {
                    currentCDS.startIdx = getLastORFPos(seq->seq.l, frameIdx);
                    currentCDS.endIdx = seq->seq.l - 1 - (frameIdx % 3);
                    currentCDS.strand = -1;
                    currentCDS.phase = 0;
            	}

				convertToAA(seq->seq.s, &currentCDS, &outputBuffer, &outputSize);

				// Line ID of CDS : Frame : Sequence Header
				fprintf(outputPtr, ">%lu:%lu:%s\n%.*s\n", currentLine, frameIdx, seq->name.s, (int) outputSize, outputBuffer);
				free(outputBuffer);

				if (!passMulti) break;
            }

            currentLine++;
    }

    // Close files
    fflush(outputPtr);
    err_gzclose(inputSeqPtr);
    err_fclose(outputPtr);
    kseq_destroy(seq);

    return 0;
}

int writeIndexDirectProtein(const char * passPrefix, const char * passProName) {
	FILE * outputPtr;

	outputPtr = fopen(passProName, "w");

	// Write index type header
	writeIndexHeader(outputPtr, 0, 0);

    // Close file
    fflush(outputPtr);
    err_fclose(outputPtr);

    return 0;
}

// This is a function use for testing - to be removed in final code
// Currently writes GC content to sequence header in 2nd position
int writeIndexTestProtein(const char * passPrefix, const char * passProName) {
    struct CDS currentCDS;
    gzFile  inputSeqPtr;
    FILE * outputPtr;
    char * outputBuffer;
    kseq_t * seq;
    int idx, gcCount;
    unsigned long outputSize, currentLine, frameIdx;

    // Prepare file handles
    inputSeqPtr = xzopen(passPrefix, "r");
    outputPtr = fopen(passProName, "w");

    // Iterate through each sequence
    seq = kseq_init(inputSeqPtr);
    currentLine = 0;

    while (kseq_read(seq) > 0) {

            gcCount = 0;
            for (idx = 0 ; idx < seq->seq.l ; idx++) {
               if (seq->seq.s[idx] == 'G' || seq->seq.s[idx] == 'C') gcCount++; 
            }
   
            gcCount = gcCount * 100 / seq->seq.l;

            for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
                    currentCDS.startIdx = (frameIdx % 3);
                    currentCDS.endIdx = seq->seq.l - 1 - 3 + (frameIdx % 3);
                    currentCDS.strand = (frameIdx < 3 ? 1 : -1);
                    currentCDS.phase = 0;

                    convertToAA(seq->seq.s, &currentCDS, &outputBuffer, &outputSize);

        			// Sequence ID : ORF Index per Sequence : Relative Frame per Sequence : Sequence Header
                    fprintf(outputPtr, ">%lu:%lu:%d:%s\n%.*s\n", currentLine, frameIdx, gcCount, seq->name.s, (int) outputSize, outputBuffer);
                    free(outputBuffer);
            }

            currentLine++;
    }

    // Close files
    fflush(outputPtr);
    err_gzclose(inputSeqPtr);
    err_fclose(outputPtr);
    kseq_destroy(seq);	

    return 0;
}


// Detects ORFs in the given nucleotide FASTA file and converts to a protein FASTA file
int writeReadsProtein(const char * passPrefix, const char * passProName, mem_opt_t * passOptions) {
	struct CDS * orfList;
	gzFile inputSeqPtr;
	FILE * outputProPtr, * outputNTPtr;
 	char * outputProBuffer, * outputNTName;
	kseq_t * seq;
	unsigned long seqIdx, outputSize, orfCount, orfIdx;

	// Check for incompatible combinations
	if ((passOptions->proteinFlag & ALIGN_FLAG_BRUTEORF) && (passOptions->indexFlag & INDEX_FLAG_MF)) {
		fprintf(stderr, "[paladin_align] Brute force ORF detection redundant to MF index, disabling...\n");
		passOptions->flag &= ~ALIGN_FLAG_BRUTEORF;
	}

	// Prepare file handles
	inputSeqPtr = xzopen(passPrefix, "r");
	outputProPtr = fopen(passProName, "w");
	outputNTPtr = 0;

	if (passOptions->proteinFlag & ALIGN_FLAG_GEN_NT) {
		outputNTName = malloc(strlen(passPrefix) + 5);
		sprintf(outputNTName, "%s.orf", passPrefix);
		outputNTPtr = fopen(outputNTName, "w");
		free(outputNTName);
	}

	// Iterate through each read
	seqIdx = 0;
	seq = kseq_init(inputSeqPtr);

	while(kseq_read(seq) >= 0) {
		// Search for ORFs
		getSequenceORF(seq->seq.s, seq->seq.l-1, passOptions, &orfList, &orfCount);

		// Write out the corresponding protein sequence for each ORF
		for (orfIdx = 0 ; orfIdx < orfCount ; orfIdx++) {
			convertToAA(seq->seq.s, orfList+orfIdx, &outputProBuffer, &outputSize);

			// Sequence ID : ORF Index per Sequence : Relative Frame per Sequence : Sequence Header
			fprintf(outputProPtr, ">%lu:%lu:%hd:%s\n%.*s\n", seqIdx, orfIdx, orfList[orfIdx].relFrame, seq->name.s, (int) outputSize, outputProBuffer);
			free(outputProBuffer);

			if (passOptions->proteinFlag & ALIGN_FLAG_GEN_NT) {
				// Sequence ID : ORF Index per Sequence : Relative Frame per Sequence : Sequence Header
				fprintf(outputNTPtr, ">%lu:%lu:%hd:%s\n%.*s\n", seqIdx, orfIdx, orfList[orfIdx].relFrame, seq->name.s, (int) (orfList[orfIdx].endIdx - orfList[orfIdx].startIdx + 1), seq->seq.s + orfList[orfIdx].startIdx);

			}
		}

		free(orfList);
		seqIdx++;
	}

	// Close files
	fflush(outputProPtr);
	err_gzclose(inputSeqPtr);
	err_fclose(outputProPtr);
	kseq_destroy(seq);

	if (passOptions->proteinFlag & ALIGN_FLAG_GEN_NT) {
		fflush(outputNTPtr);
		err_fclose(outputNTPtr);
	}

	return 0;
}
