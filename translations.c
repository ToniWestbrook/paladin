#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "translations.h"

// Codon encoded as 6 bit value, MSB as left-most nucleotide in 3-mer
unsigned char codon_aa_hash[][64] = {
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', '*', 'S', '*', 'S', 'M', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'M', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'T', 'T', 'T', 'T', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'S', 'S', 'S', 'S', 'M', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', 'Q', 'Y', 'Q', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*'}, 
	{'*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*'}, 
	{'N', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'S', 'S', 'S', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'C', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'S', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'G', 'S', 'G', 'S', 'M', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'N', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'S', 'S', 'S', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', 'Y', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', 'Q', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', 'L', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*'}, 
	{'*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*'}, 
	{'*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*'}, 
	{'*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*'}, 
	{'N', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'S', 'S', 'S', 'S', 'M', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', 'L', 'Y', '*', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', '*', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'S', 'S', 'K', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'W', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}, 
	{'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', 'G', 'C', 'W', 'C', 'L', 'F', 'L', 'F'}
};

const char * translationNames[] = {"Standard", "Vertebrate Mitochondrial", "Yeast Mitochondrial", "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma", "Invertebrate Mitochondrial", "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear", "", "", "Echinoderm Mitochondrial; Flatworm Mitochondrial", "Euplotid Nuclear", "Bacterial, Archaeal and Plant Plastid", "Alternative Yeast Nuclear", "Ascidian Mitochondrial", "Alternative Flatworm Mitochondrial", "Blepharisma Macronuclear", "Chlorophycean Mitochondrial", "", "", "", "", "Trematode Mitochondrial", "Scenedesmus obliquus Mitochondrial", "Thraustochytrium Mitochondrial", "Pterobranchia Mitochondrial", "Candidate Division SR1 and Gracilibacteria"};

// Convert command line comma delimited translation string into int array
int * convertTransArgs(const char * passArguments) {
    char * workArguments, * parseArgument;
    int * retTranslations;
    int argIdx;
    int transCount, valid;

    // Check for help mode
    if (passArguments[0] == '?') {
        renderTranslations();
        exit(0);
    }

    // Copy arguments over to temporary area
    workArguments = calloc(strlen(passArguments)+1, 1);
    memcpy(workArguments, passArguments, strlen(passArguments)+1);

    // Count commas to calculate array size, and replace with NULLs for parsing
    for (argIdx = 0, transCount = 0 ; workArguments[argIdx] ; argIdx++) {
        if (workArguments[argIdx] == ',') {
            workArguments[argIdx] = 0;
            transCount++;
        }
    }

    transCount++;
    retTranslations = calloc(transCount+1, sizeof(int));

    // Parse arguments
    parseArgument = workArguments;
    valid = 1;

    for (argIdx = 0 ; argIdx < transCount ; argIdx++) {
        retTranslations[argIdx] = atoi(parseArgument); 
        parseArgument += strlen(parseArgument) + 1;

        // Check for valid translation ID
        if (retTranslations[argIdx] < 1 || retTranslations[argIdx] > TRANSLATION_COUNT) {
            valid = 0; break;
        }

        if (translationNames[retTranslations[argIdx] - 1][0] == 0) {
            valid = 0; break;
        }
    }

    // If invalid, only translate for standard genetic code
    if (!valid) {
        retTranslations[0] = 1;
        retTranslations[1] = 0;
    }

    return retTranslations;
}

// Print possible translations
void renderTranslations() {
    char * transName;
    int transIdx;

    printf("\nThe following genetic code translation tables are available to PALADIN.\n"
           "Multiple codes may be specified by a comma-delimited list with no spaces.\n"
           "If no translation codes are specified, code 1 is used.\n"
           "Details may be found at http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi\n\n"); 

    for (transIdx = 0 ; transIdx < TRANSLATION_COUNT ; transIdx++) {
        if ((transName = (char *) translationNames[transIdx])[0] != 0) {
            printf("%d. %s\n", transIdx + 1, transName); 
        }
    }
}
