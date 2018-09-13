#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include "protein.h"
#include "translations.h"
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

/*unsigned char codon_aa_hash[1][64] = {{
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
}};
*/
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
int convertToAA(char * passSequence, CDS * passCDS, int passTrans, char ** retSequence, unsigned long * retSize) {
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

		if (currentCodon == 0xFF) {
            // If ambiguous nucleotides are present, assign ambiguous AA
			(*retSequence)[aaIdx] = 'X';
		}
        else {
            // Hash to amino acid IUPAC
            (*retSequence)[aaIdx] = codon_aa_hash[passTrans][currentCodon];
        }
	}

	return 0;
}

// Calculate last aligned position in sequence
long getLastAlignedPos(long passLength, int passFrame) {
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

long getLastAlignedOrfPos(long passLength, int passFrame, mem_opt_t * passOptions) {
	long retValue;

	if (passOptions->min_orf_percent) {
		// Minimum ORF length specified as a read length percentage
		retValue = getLastAlignedPos(passOptions->min_orf_percent * passLength, passFrame);
	}
	else {
		// Minimum ORF length specified as constant value, check for adjustment
		if ((passOptions->proteinFlag & ALIGN_FLAG_ADJUST_ORF) && (passLength < passOptions->min_orf_len)) {
			retValue = getLastAlignedPos(passLength, passFrame);
		}
		else {
			retValue = getLastAlignedPos(passOptions->min_orf_len, passFrame);
		}
	}


	return retValue;
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
void compileORFHistory(long * passHistory[2][6], long passHistorySize[6], CDS * * retCDS, unsigned long * retCount) {
	unsigned long srcFrameIdx, srcHistoryIdx;
	unsigned long srcStart, srcEnd, totalSize;
	int relStart, validEntry;

	// Allocate potential CDS records per entry in history array
	totalSize = 0;
	for (srcFrameIdx = 0 ; srcFrameIdx < 6 ; srcFrameIdx++) totalSize += passHistorySize[srcFrameIdx];

	*retCDS = calloc(totalSize, sizeof(CDS));
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
int getSequenceORF(char * passSequence, unsigned long passLength, int passTrans, mem_opt_t * passOptions, CDS * * retCDS, unsigned long * retCount) {
	int frameIdx, addIdx, strandDir, relFrame;
	long seqIdx, absIdx, lastStart;
	long endSeqPos, endOrfPos;
	unsigned char currentCodon;
	long * history[2][6], historySize[6];

	// Initialize values
	*retCDS = NULL;
	*retCount = 0;
	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		history[0][frameIdx] = NULL;
		history[1][frameIdx] = NULL;
		historySize[frameIdx] = 0;
	}

	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		strandDir = (frameIdx < 3 ? 1 : -1);
		relFrame = frameIdx % 3;
		lastStart = relFrame;
		endSeqPos = getLastAlignedPos(passLength, relFrame);
		endOrfPos = getLastAlignedOrfPos(passLength, relFrame, passOptions);

		// Adjust min ORF length if requested
		if ((passOptions->proteinFlag & ALIGN_FLAG_ADJUST_ORF) && (passLength < passOptions->min_orf_len)) {
			endOrfPos = getLastAlignedPos(passLength, relFrame);
		}

		// Iterate through all codons in this read frame
		for (seqIdx = relFrame ; seqIdx + 2 <= endSeqPos ; seqIdx += 3) {
			// Calculate absolute index from relative
			absIdx = seqIdx;
			if (strandDir == -1) absIdx = absIdx * -1 + passLength - 1;

			// Encode codon
			currentCodon = encodeCodon(passSequence + absIdx, strandDir);

			// Check for stop codon or EOS
			if ((codon_aa_hash[passTrans][currentCodon] == '*') || (seqIdx + 2 == endSeqPos)) {
				if (seqIdx + relFrame + 2 - lastStart >= endOrfPos) {
					// Translate current (non-brute) or all frames (brute)
					for (addIdx = 0 ; addIdx < 6 ; addIdx++) {
						if ((passOptions->proteinFlag & ALIGN_FLAG_BRUTE_ORF) || (addIdx == frameIdx)) {
							strandDir = (addIdx < 3 ? 1 : -1);
							relFrame = addIdx % 3;

							addORFHistory(history, historySize, addIdx);
							history[0][addIdx][historySize[addIdx] - 1] = (strandDir == 1) ? relFrame : passLength - 1 - relFrame;
							history[1][addIdx][historySize[addIdx] - 1] = getLastAlignedPos(passLength, addIdx);
						}
					}

					break;
				}

				lastStart = seqIdx + 3;
			}
		}

		// If we have brute-force match, end processing
		if (frameIdx < 5) {
			if (historySize[frameIdx + 1]) break;
		}
	}

	// Convert history array into CDS entries
	compileORFHistory(history, historySize, retCDS, retCount);

	for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
		free(history[0][frameIdx]);
		free(history[1][frameIdx]);
	}

	return 0;
}

// Iterate through GFF annotation file, return next available CDS entry
int getNextCDS(FILE * passFile, CDS * retCDS, unsigned long * retLine) {
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
int writeIndexProtein(const char * passPrefix, const char * passProName, const char * passAnnName, IndexHeader passHeader) {
	CDS currentCDS;
	gzFile inputSeqPtr;
	FILE * inputAnnPtr, * outputPtr;
	char * outputBuffer;
	kseq_t * seq;
	unsigned long outputSize, currentLine, frameIdx;

	// Prepare file handles
	inputSeqPtr = xzopen(passPrefix, "r");
	inputAnnPtr = err_xopen_core(__func__, passAnnName, "r");
	outputPtr = err_xopen_core(__func__, passProName, "w");

	// Write index type header
	passHeader.nucleotide = 1;
	writeIndexHeader(outputPtr, passHeader);

	// Read in 1st sequence data
	seq = kseq_init(inputSeqPtr);
	kseq_read(seq);

	// Iterate through each CDS sequence in the annotation file
	currentLine = 0;
	while (getNextCDS(inputAnnPtr, &currentCDS, &currentLine)) {
		for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
			convertToAA(seq->seq.s, &currentCDS, 0, &outputBuffer, &outputSize);

			// Line ID of CDS : Frame : Sequence Header
			fprintf(outputPtr, ">%lu:%lu:%s\n%.*s\n", currentLine, frameIdx, currentCDS.description, (int) outputSize, outputBuffer);
			free(outputBuffer);

			currentCDS.startIdx++;
			currentCDS.endIdx++;
			if (frameIdx == 2) currentCDS.strand *= -1;

			if (!passHeader.multiFrame) break;
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

int writeIndexCodingProtein(const char * passPrefix, const char * passProName, IndexHeader passHeader) {
	CDS currentCDS;
	gzFile  inputSeqPtr;
	FILE * outputPtr;
	char * outputBuffer;
	kseq_t * seq;
	unsigned long outputSize, currentLine, frameIdx;

	// Prepare file handles
	inputSeqPtr = xzopen(passPrefix, "r");
	outputPtr = err_xopen_core(__func__, passProName, "w");

	// Write index type header
	passHeader.nucleotide = 1;
	writeIndexHeader(outputPtr, passHeader);

	// Iterate through each sequence
	seq = kseq_init(inputSeqPtr);
	currentLine = 0;

	while (kseq_read(seq) > 0) {
		for (frameIdx = 0 ; frameIdx < 6 ; frameIdx++) {
			if (frameIdx < 3) {
				currentCDS.startIdx = (frameIdx % 3);
				currentCDS.endIdx = getLastAlignedPos(seq->seq.l, frameIdx);
				currentCDS.strand = 1;
				currentCDS.phase = 0;
			}
			else {
				currentCDS.startIdx = getLastAlignedPos(seq->seq.l, frameIdx);
				currentCDS.endIdx = seq->seq.l - 1 - (frameIdx % 3);
				currentCDS.strand = -1;
				currentCDS.phase = 0;
			}

			convertToAA(seq->seq.s, &currentCDS, 0, &outputBuffer, &outputSize);

			// Line ID of CDS : Frame : Sequence Header
			fprintf(outputPtr, ">%lu:%lu:%s\n%.*s\n", currentLine, frameIdx, seq->name.s, (int) outputSize, outputBuffer);
			free(outputBuffer);

			if (!passHeader.multiFrame) break;
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

int writeIndexDirectProtein(const char * passPrefix, const char * passProName, IndexHeader passHeader) {
	FILE * outputPtr;

	outputPtr = err_xopen_core(__func__, passProName, "w");

	// Write index type header
	passHeader.nucleotide = 0;
	writeIndexHeader(outputPtr, passHeader);

	// Close file
	fflush(outputPtr);
	err_fclose(outputPtr);

	return 0;
}

// This is a function use for testing - to be removed in final code
// Currently writes GC content to sequence header in 2nd position
int writeIndexTestProtein(const char * passPrefix, const char * passProName) {
	CDS currentCDS;
	gzFile  inputSeqPtr;
	FILE * outputPtr;
	char * outputBuffer;
	kseq_t * seq;
	int idx, gcCount;
	unsigned long outputSize, currentLine, frameIdx;

	// Prepare file handles
	inputSeqPtr = xzopen(passPrefix, "r");
	outputPtr = err_xopen_core(__func__, passProName, "w");

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

			convertToAA(seq->seq.s, &currentCDS, 0, &outputBuffer, &outputSize);

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
	CDS * orfList;
	gzFile inputSeqPtr;
	FILE * outputProPtr, * outputNTPtr;
	char * outputProBuffer, * outputNTName;
	kseq_t * seq;
	unsigned long seqIdx, transIdx, transTable, outputSize, orfCount, orfTotal, orfIdx;

	// Check for incompatible combinations
	if ((passOptions->proteinFlag & ALIGN_FLAG_BRUTE_ORF) && passOptions->indexInfo.multiFrame) {
		logMessage(__func__, LOG_LEVEL_WARNING, "Brute force ORF detection redundant to MF index, disabling...\n");
		passOptions->flag &= ~ALIGN_FLAG_BRUTE_ORF;
	}

	// Prepare file handles
	inputSeqPtr = xzopen(passPrefix, "r");
	outputProPtr = err_xopen_core(__func__, passProName, "w");
	outputNTPtr = 0;

	if (passOptions->proteinFlag & ALIGN_FLAG_GEN_NT) {
		outputNTName = malloc(strlen(passPrefix) + 5);
		sprintf(outputNTName, "%s.orf", passPrefix);
		outputNTPtr = err_xopen_core(__func__, outputNTName, "w");
		free(outputNTName);
	}

	logMessage(__func__, LOG_LEVEL_MESSAGE, "Detecting open reading frames...\n");

	// Iterate through each read
	orfTotal = 0;
	seqIdx = 0;
	seq = kseq_init(inputSeqPtr);

	while(kseq_read(seq) >= 0) {
        transIdx = 0;

        while((passOptions->translations)[transIdx++] != 0) {
            // Lookup current translation table
            transTable = (passOptions->translations)[transIdx - 1] - 1;

            // Search for ORFs
            getSequenceORF(seq->seq.s, seq->seq.l, transTable, passOptions, &orfList, &orfCount);

            // Write out the corresponding protein sequence for each ORF
            for (orfIdx = 0 ; orfIdx < orfCount ; orfIdx++) {
                convertToAA(seq->seq.s, orfList+orfIdx, transTable, &outputProBuffer, &outputSize);

                // Sequence ID : ORF Index per Sequence : Relative Frame per Sequence : Sequence Header
                fprintf(outputProPtr, ">%lu:%lu:%lu:%s\n%.*s\n", 
                        seqIdx, orfIdx, transTable * 6 + orfList[orfIdx].relFrame, seq->name.s, (int) outputSize, outputProBuffer);
                free(outputProBuffer);

                if (passOptions->proteinFlag & ALIGN_FLAG_GEN_NT) {
                    // Sequence ID : ORF Index per Sequence : Relative Frame per Sequence : Sequence Header
                    fprintf(outputNTPtr, ">%lu:%lu:%lu:%s\n%.*s\n", 
                            seqIdx, orfIdx, transTable * 6 + orfList[orfIdx].relFrame, seq->name.s, (int) (orfList[orfIdx].endIdx - orfList[orfIdx].startIdx + 1), seq->seq.s + orfList[orfIdx].startIdx);

                }
            }

            free(orfList);
            orfTotal += orfCount;
        }
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

	if (passOptions->proteinFlag & ALIGN_FLAG_BRUTE_ORF) orfTotal /= 6;
	logMessage(__func__, LOG_LEVEL_MESSAGE, "Detected and translated %d open reading frames in %d sequences\n", orfTotal, seqIdx);

	return 0;
}
