# PALADIN
Protein Alignment and Detection Interface

 PALADIN is a protein sequence alignment tool based on the BWA source.  Like BWA, it aligns
 sequences via read-mapping using BWT.  PALADIN, however, offers the novel approach
 of aligning in the protein space.  During the index phase, it processes the reference genome's
 nucleotide sequences and GTF/GFF annotation containing CDS entries, first
 converting these transcripts into the corresponding protein sequences, then creating the BWT
 and suffix array from these proteins.  During the alignment phase, it attempts to find ORFs in 
 the read sequences, then converts these to protein sequences, and aligns to the reference 
 protein sequences. 

 PALADIN currently only supports single-end reads, and BWA-MEM based alignment.  It makes 
 use of many BWA parameters and is therefore compatible with many of its command line arguments.
