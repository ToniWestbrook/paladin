/* Contact: Toni Westbrook <anthonyw@wildcats.unh.edu> 
*
* PALADIN (Protein Alignment and Detection Interface)
*
* PALADIN is a protein sequence alignment tool based on the BWA source.  Like BWA, it aligns
* sequences via read-mapping using BWT.  PALADIN, however, offers the novel approach
* of aligning in the protein space.  During the index phase, it processes the reference genome's
* nucleotide sequences and GTF/GFF annotation containing CDS entries, first
* converting these transcripts into the corresponding protein sequences, then creating the BWT
* and suffix array from these proteins.  During the alignment phase, it attempts to find ORFs in 
* the read sequences, then converts these to protein sequences, and aligns to the reference 
* protein sequences. 
*
* PALADIN currently only supports single-end reads, and BWA-MEM based alignment.  It makes 
* use of many BWA parameters and is therefore compatible with many of its command line arguments.
*
*
* PALADIN IS CURRENTLY PRE-ALPHA AND HAS NOT BEEN FULLY TESTED.  USE AT YOUR OWN RISK.  
* 
*
* For information regarding BWA, please contact its author, Heng Li <lh3@sanger.ac.uk> */ 

#include <stdio.h>
#include <string.h>
#include "kstring.h"
#include "utils.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.0"
#endif

int bwa_fa2pac(int argc, char *argv[]);
int bwa_pac2bwt(int argc, char *argv[]);
int bwa_bwtupdate(int argc, char *argv[]);
int bwa_bwt2sa(int argc, char *argv[]);
int bwa_index(int argc, char *argv[]);

int main_mem(int argc, char *argv[]);
int main_shm(int argc, char *argv[]);

int main_pemerge(int argc, char *argv[]);

int writeProtein(const char * passPrefix, const char * passProName, const char * passAnnName);

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: PALADIN (Protein Alignment and Detection Interface)\n");

	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Toni Westbrook <anthonyw@wildcats.unh.edu>\n");
	fprintf(stderr, "Based on: BWA by Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   paladin <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in FASTA format with GTF/GFF annotations\n");
	fprintf(stderr, "         mem           BWA-MEM algorithm\n");
	fprintf(stderr, "         pemerge       merge overlapping paired ends (EXPERIMENTAL)\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         shm           manage indices in shared memory\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");
	fprintf(stderr, "\n");
	fprintf(stderr,
"Note: To use PALADIN, you must first index the genome with `paladin index'.\n"
"      Alignment may then be performed using the PALADIN's BWA-MEM based algorithm: `mem'\n\n");
	return 1;
}

int main(int argc, char *argv[])
{
	extern char *bwa_pg;
	int i, ret;
	double t_real;
	kstring_t pg = {0,0,0};
	t_real = realtime();

	ksprintf(&pg, "@PG\tID:PALADIN\tPN:PALADIN\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);

	// Process command argument
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "fa2pac") == 0) ret = bwa_fa2pac(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwt") == 0) ret = bwa_pac2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtupdate") == 0) ret = bwa_bwtupdate(argc-1, argv+1);
	else if (strcmp(argv[1], "bwt2sa") == 0) ret = bwa_bwt2sa(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) ret = bwa_index(argc-1, argv+1);
	else if (strcmp(argv[1], "mem") == 0) ret = main_mem(argc-1, argv+1);
	else if (strcmp(argv[1], "shm") == 0) ret = main_shm(argc-1, argv+1);
	else if (strcmp(argv[1], "pemerge") == 0) ret = main_pemerge(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	err_fflush(stdout);
	err_fclose(stdout);

	// If execution was successful, print time and CPU statistics
	if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}
	free(bwa_pg);
	return ret;
}
