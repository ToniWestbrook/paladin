/*
   The MIT License

   Copyright (c) 2015 by Anthony Westbrook, University of New Hampshire <anthony.westbrook@unh.edu>
   Copyright (c) 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
   PALADIN (Protein Alignment and Detection Interface)

   PALADIN is a protein sequence alignment tool designed for the accurate
   functional characterization of metagenomes.

   PALADIN is based on BWA, and aligns sequences via read-mapping using
   BWT. PALADIN, however, offers the novel approach of aligning in the
   protein space. During the index phase, it processes the reference genome's
   nucleotide sequences and GTF/GFF annotation containing CDS entries, first
   converting these transcripts into the corresponding protein sequences, then
   creating the BWT and suffix array from these proteins. The process of
   translation is skiped when providing a protein reference file (e.g., UniProt)
   for mapping. During the alignment phase, it attempts to find ORFs in the
   read sequences, then converts these to protein sequences, and aligns to the
   reference protein sequences.

   PALADIN currently only supports single-end reads (or reads merged with FLASH,
   PEAR, abyss-mergepairs), and BWA-MEM based alignment. It makes use of many
   BWA parameters and is therefore compatible with many of its command line
   arguments.

   PALADIN may output a standard SAM file, or a text file containing a
   UniProt-generated functional profile. This text file may be used for all
   downstream characterizations.

   Contact: Toni Westbrook <anthony.westbrook@unh.edu>
   For information regarding BWA, contact Heng Li <lh3@sanger.ac.uk>
*/

#include <stdio.h>
#include <string.h>
#include "main.h"
#include "align.h"
#include "bwtindex.h"
#include "kstring.h"
#include "utils.h"

int main(int argc, char *argv[]) {
	extern char *bwa_pg;
	int i, ret;
	double t_real;
	kstring_t pg = {0,0,0};
	t_real = realtime();

	ksprintf(&pg, "@PG\tID:PALADIN\tPN:PALADIN\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);

	// Process command argument
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;

	if (argc < 2) return renderMainUsage();
	if (strcmp(argv[1], "index") == 0) ret = command_index(argc-1, argv+1);
	else if (strcmp(argv[1], "prepare") == 0) ret = command_prepare(argc-1, argv+1);
	else if (strcmp(argv[1], "align") == 0) ret = command_align(argc-1, argv+1);
	else if (strcmp(argv[1], "fa2pac") == 0) ret = bwa_fa2pac(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwt") == 0) ret = command_pac2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtupdate") == 0) ret = command_bwtupdate(argc-1, argv+1);
	else if (strcmp(argv[1], "bwt2sa") == 0) ret = command_bwt2sa(argc-1, argv+1);
	else if (strcmp(argv[1], "shm") == 0) ret = main_shm(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) return renderVersion();
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

int renderMainUsage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: PALADIN (Protein Alignment and Detection Interface)\n\n");
	fprintf(stderr, "Usage:   paladin <command> [options]\n\n");
	fprintf(stderr, "Command: index         index NT or AA sequences in FASTA format\n");
	fprintf(stderr, "         prepare       download and index protein reference\n");
	fprintf(stderr, "         align         align single end read sequences\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         shm           manage indices in shared memory\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n\n");
	fprintf(stderr, "         version       version and contact information\n");
	fprintf(stderr, "\n");
	fprintf(stderr,
"Note: To use PALADIN, first index the reference using the 'index' or 'prepare' command.\n"
"      Then align your reads using the 'align' command.\n\n");

	return 1;
}

int renderVersion() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: PALADIN (Protein Alignment and Detection Interface)\n");

	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Toni Westbrook (UNH) <anthony.westbrook@unh.edu>\n");
	fprintf(stderr, "Based on: BWA by Heng Li <lh3@sanger.ac.uk>\n\n");

	return 1;
}
