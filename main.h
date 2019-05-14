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


#ifndef MAIN_H_
#define MAIN_H_

#define STR_CONVERT(x) #x
#define STR(x) STR_CONVERT(x)

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION STR(PACKAGE_VERSION_MAJOR) "." STR(PACKAGE_VERSION_MINOR) "." STR(PACKAGE_VERSION_REV)
#define PACKAGE_VERSION_MAJOR 1
#define PACKAGE_VERSION_MINOR 4
#define PACKAGE_VERSION_REV 6
#endif

// Render usage and version details
int renderMainUsage();
int renderVersion();

// CLEAN
int bwa_fa2pac(int argc, char *argv[]);
int main_shm(int argc, char *argv[]);
int main_pemerge(int argc, char *argv[]);

#endif /* MAIN_H_ */
