/* Contact: Toni Westbrook <anthonyw@wildcats.unh.edu> */

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <unistd.h>
#include "protein.h"
#include "bntseq.h"
#include "utils.h"
#include "bwt.h"

#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_MAP_INIT_STR(str, int)

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

// Construct 6-bit codon from 3 ASCII characters
unsigned char encodeCodon(unsigned char * passSequence, int passStrand) {
	unsigned char retCodon;
	int relIdx, absIdx;

	retCodon = 0;

	for (relIdx = 0 ; relIdx < 3 ; relIdx++) {
		// Calculate absolute index from relative
		absIdx = relIdx;
		if (passStrand < 0) absIdx = absIdx * -1;
 
		// If this codon has any ambiguous nucleotides, signal with 0xFF
		if (nst_nt4_table[passSequence[absIdx]] > 3) {
			return 0xFF;
		}

		retCodon = (retCodon << 2) | nst_nt4_table[passSequence[absIdx]];
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
		
		// If ambiguous nucleotides are present, assign ambiguous AA and exit
		if (currentCodon == 0xFF) {
			(*retSequence)[aaIdx] = 'X';
			break;
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
	unsigned long srcFrameIdx, srcHistoryIdx, dstFrameIdx, dstHistoryIdx, totalSize;
	unsigned long srcStart, srcEnd, dstStart, dstEnd, srcLen, dstLen;
	int validEntry;

	// Start by allocating potential CDS records per entry in history array
	totalSize = 0;
	for (srcFrameIdx = 0 ; srcFrameIdx < 6 ; srcFrameIdx++) totalSize += passHistorySize[srcFrameIdx];
	*retCDS = calloc(totalSize, sizeof(struct CDS));
	*retCount = 0;

	// Iterate through each entry in each reading frame
	for (srcFrameIdx = 0 ; srcFrameIdx < 6 ; srcFrameIdx++) {
		for (srcHistoryIdx = 0 ; srcHistoryIdx < passHistorySize[srcFrameIdx] ; srcHistoryIdx++) {
            // Swap start/stop for reverse direction
			if (srcFrameIdx < 3) {
				srcStart = passHistory[0][srcFrameIdx][srcHistoryIdx];
				srcEnd = passHistory[1][srcFrameIdx][srcHistoryIdx];
			}
			else {
				srcStart = passHistory[1][srcFrameIdx][srcHistoryIdx];
				srcEnd = passHistory[0][srcFrameIdx][srcHistoryIdx];
			}

			srcLen = srcEnd - srcStart + 1;
			validEntry = 1;

			// Filter minimum ORF size (to be set from command line args)
			if (srcLen < 180) continue;
			
			// Search for larger overlapping ORFs
			for (dstFrameIdx = 0 ; dstFrameIdx < 6 ; dstFrameIdx++) {
				// No need to compare same frame
				if (srcFrameIdx == dstFrameIdx) continue;

				// Check for larger, overlapping frames
				for (dstHistoryIdx = 0 ; dstHistoryIdx < passHistorySize[dstFrameIdx] ; dstHistoryIdx++) {
					if (dstFrameIdx < 3) {
						dstStart = passHistory[0][dstFrameIdx][dstHistoryIdx];
						dstEnd = passHistory[1][dstFrameIdx][dstHistoryIdx];
					}
					else {
						dstStart = passHistory[1][dstFrameIdx][dstHistoryIdx];
						dstEnd = passHistory[0][dstFrameIdx][dstHistoryIdx];
					}

					dstLen = abs(dstEnd - dstStart) + 1;

					// Compare potential frames that are larger 
					if (dstLen  > srcLen) {
						// src.start <= dst.start <= src.end
						if ((srcStart + 0 <= dstStart) && (dstStart <= srcEnd - 0)) validEntry = 0;

						// dst.start <= src.start <= dst.end
						if ((dstStart + 0  <= srcStart) && (srcStart <= dstEnd - 0)) validEntry = 0;

						if (!validEntry) break;
					}
				}

				if (!validEntry) break;
			}

			// Add CDS record if valid
			if (validEntry) {
				(*retCDS)[*retCount].startIdx = srcStart;
				(*retCDS)[*retCount].endIdx = srcEnd;
				(*retCDS)[*retCount].strand = (srcFrameIdx < 3) ? 1 : -1;
				
				(*retCount)++;
			}
		}
	}
}

// Scan nucleotide sequence for all recognized ORFs, return as CDS array
int getSequenceORF(char * passSequence, unsigned long passLength, struct CDS * * retCDS, unsigned long * retCount) {
	unsigned long strandIdx, frameIdx, seqIdx, historyIdx;
	long * history[2][6], historySize[6], absIdx;
	unsigned char currentCodon;
	int frameStart;

	for (strandIdx = 0 ; strandIdx < 2 ; strandIdx++) {
		for (frameIdx = 0 ; frameIdx < 3 ; frameIdx++) {
			// Initialize history (start, end)
			historyIdx = strandIdx * 3 + frameIdx;
			historySize[historyIdx] = 0;
			frameStart = 0;
			
			// Iterate through all codons in this read frame
			for (seqIdx = frameIdx ; seqIdx + 3 < passLength ; seqIdx += 3) {
				// Calculate absolute index from relative
				absIdx = seqIdx;
				if (strandIdx == 1) absIdx = absIdx * -1 + passLength - 1;

				// Encode codon
				currentCodon = encodeCodon(passSequence + absIdx, (strandIdx == 0) ? 1 : -1);

				// Scan for start codon (ATG - 0x0E)
				if (currentCodon == 0x0E) {
					if (!frameStart) {
						addORFHistory(history, historySize, historyIdx);

						history[0][historyIdx][historySize[historyIdx] - 1] = absIdx;
						history[1][historyIdx][historySize[historyIdx] - 1] = getLastORFPos(passLength, strandIdx * 3 + frameIdx) ;
						frameStart = 1;
					}
				}
 
				// Scan for stop codon (TAA - 0x30, TAG - 0x32, TGA - 0x38)
				if ((currentCodon == 0x30) || (currentCodon == 0x32) || (currentCodon == 0x38)) {
					// Check if our ORF started before the sequence
					if (!frameStart && (historySize[historyIdx] == 0)) {
						addORFHistory(history, historySize, historyIdx);
						history[0][historyIdx][historySize[historyIdx] - 1] = (strandIdx == 0) ? frameIdx : passLength - 1 - frameIdx;
						frameStart = 1;
					}

					// Ignore stop codons outside of ORFs
					if (frameStart) {
						history[1][historyIdx][historySize[historyIdx] - 1] = absIdx + ((strandIdx == 0) ? 2 : -2);
						frameStart = 0;
					}
				}
			}

			// Check if entire frame was open
			if (historySize[historyIdx] == 0) {
				addORFHistory(history, historySize, historyIdx);
				history[0][historyIdx][0] = (strandIdx == 0) ? frameIdx : (passLength - 1 - frameIdx);
				history[1][historyIdx][0] = getLastORFPos(passLength, strandIdx * 3 + frameIdx);
			}
		}
	}

	// Filter, remove overlaps, create CDS array
	compileORFHistory(history, historySize, retCDS, retCount);

	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) free(history[frameIdx / 3][frameIdx % 3]);

	return 0;
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

			// Quick code to insert description, fix this up
			sprintf(retCDS->description, "%s", field[8]);
			int charIdx;
			for (charIdx = 0 ; charIdx < strlen(retCDS->description) ; charIdx++) {
				if (retCDS->description[charIdx] == ' ') retCDS->description[charIdx] = '_';
			}

			return 1;
		}
	}

	return 0;
}

// Convert the given reference nucleotide FASTA and GFF file to a protein FASTA file
int writeIndexProtein(const char * passPrefix, const char * passProName, const char * passAnnName) {
	struct CDS currentCDS;
	gzFile  inputSeqPtr;
	FILE * inputAnnPtr, * outputPtr;
 	char * outputBuffer;
	kseq_t * seq;
	unsigned long outputSize, currentLine;

	// Prepare file handles
	inputSeqPtr = xzopen(passPrefix, "r");
	inputAnnPtr = fopen(passAnnName, "r");
	outputPtr = fopen(passProName, "w");

	// Read in 1st sequence data
	seq = kseq_init(inputSeqPtr);
	kseq_read(seq);

	// Iterate through each CDS sequence in the annotation file
	currentLine = 0;
	while (getNextCDS(inputAnnPtr, &currentCDS, &currentLine)) {
		convertToAA(seq->seq.s, &currentCDS, &outputBuffer, &outputSize);
		fprintf(outputPtr, ">%d\n%.*s\n", currentLine,  outputSize, outputBuffer);
		free(outputBuffer);
	}

	// Close files
	fflush(outputPtr);
	err_gzclose(inputSeqPtr);
	err_fclose(inputAnnPtr);
	err_fclose(outputPtr);

	return 0;
}

// Convert the given reference nucleotide FASTA and GFF file to a protein FASTA file (for all 6 frames) - experimental
int writeIndexMultiProtein(const char * passPrefix, const char * passProName, const char * passAnnName) {
        struct CDS currentCDS;
        gzFile  inputSeqPtr;
        FILE * inputAnnPtr, * outputPtr;
        char * outputBuffer;
        kseq_t * seq;
        unsigned long outputSize, currentLine;

        // Prepare file handles
        inputSeqPtr = xzopen(passPrefix, "r");
        inputAnnPtr = fopen(passAnnName, "r");
        outputPtr = fopen(passProName, "w");

        // Read in 1st sequence data
        seq = kseq_init(inputSeqPtr);
        kseq_read(seq);

        // Iterate through each CDS sequence in the annotation file
        currentLine = 0;
        int test;
        while (getNextCDS(inputAnnPtr, &currentCDS, &currentLine)) {
                for (test = 0 ; test < 6 ; test++) {
                        currentCDS.startIdx++;
			currentCDS.endIdx++;
                        if (test == 3) currentCDS.strand *= -1;
			convertToAA(seq->seq.s, &currentCDS, &outputBuffer, &outputSize);
			fprintf(outputPtr, ">%d\n%.*s\n", currentLine,  outputSize, outputBuffer);
			free(outputBuffer);
		}
        }

        // Close files
        fflush(outputPtr);
        err_gzclose(inputSeqPtr);
        err_fclose(inputAnnPtr);
        err_fclose(outputPtr);

        return 0;
}

// Detects ORFs in the given nucleotide FASTA file and converts to a protein FASTA file
int writeReadsProtein(const char * passPrefix, const char * passProName) {
	struct CDS * orfList;
	gzFile inputSeqPtr;
	FILE * outputPtr;
 	char * outputBuffer;
	kseq_t * seq;
	unsigned long seqIdx, outputSize, orfCount, orfIdx;
	char testHeader[4096];

	// Prepare file handles
	inputSeqPtr = xzopen(passPrefix, "r");
	outputPtr = fopen(passProName, "w");

	// Iterate through each read
	seqIdx = 0;
	seq = kseq_init(inputSeqPtr);

	while(kseq_read(seq) >= 0) {
		// Search for ORFs
		getSequenceORF(seq->seq.s, seq->seq.l, &orfList, &orfCount);
		//if (orfCount > 0) testReadHeader(seq->name.s, testHeader);
		// Write out the corresponding protein sequence for each ORF
		for (orfIdx = 0 ; orfIdx < orfCount ; orfIdx++) {
			convertToAA(seq->seq.s, orfList+orfIdx, &outputBuffer, &outputSize);
			fprintf(outputPtr, ">%d:%d\n%.*s\n", seqIdx, orfIdx, outputSize, outputBuffer);
			free(outputBuffer);
		}

		seqIdx++;
	}

	// Close files
	fflush(outputPtr);
	err_gzclose(inputSeqPtr);
	err_fclose(outputPtr);

	return 0;
}