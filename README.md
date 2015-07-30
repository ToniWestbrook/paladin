# PALADIN

**P**rotein **AL**ignment **A**nd **D**etection **IN**terface

PALADIN is a protein sequence alignment tool designed for the accurate functional characterization of metagenomes.

PALADIN is based on BWA, and aligns sequences via read-mapping using BWT. PALADIN, however, offers the novel approach of aligning in the protein space.  During the index phase, it processes the reference genome's nucleotide sequences and GTF/GFF annotation containing CDS entries, first converting these transcripts into the corresponding protein sequences, then creating the BWT and suffix array from these proteins. The process of translatation is skiped when providing a protein reference file (e.g., UniProt) for mapping. During the alignment phase, it attempts to find ORFs in the read sequences, then converts these to protein sequences, and aligns to the reference protein sequences. 

PALADIN currently only supports single-end reads (or reads merged with FLASH, PEAR, abyss-mergepairs), and BWA-MEM based alignment. It makes use of many BWA parameters and is therefore compatible with many of its command line arguments.

PALADIN may output a standard SAM file, or a text file containing a UniProt-generated functional profile. This text file may be used for all downstream characterizations. 


INSTALLATION
--
**Dependencies**

- From a fresh install of Ubuntu, you will need to install `build-essential libcurl4-openssl-dev git make gcc zlib1g-dev`. This should be available on Ubuntu 14.04 using `sudo apt-get install build-essential libcurl4-openssl-dev git make gcc zlib1g-dev`
- PALADIN compiles by default on OSX 10.10.x

```
git clone https://github.com/twestbrookunh/paladin.git
cd paladin/
make
PATH=$PATH:$(pwd)
```

SAMPLE COMMANDS
--

Download and prepare UniProt index files.
```
paladin prepare -r0 
```

Index UniProt (or another protein) fasta, if not using the automated `prepare` command
```
paladin index -f -r2 uniprot_sprot.fasta.gz
```
Align a set of reads using 4 theads. Send the full UniProt report to paladin_uniprot_report.txt.
```
paladin align -t 4 -u 2 index input.fastq.gz > paladin_uniprot_report.txt
```
Align a set of reads using 4 theads. Produce a bam file.
```
paladin align -t 4 index input.fastq.gz | samtools view -Sb - > test.bam
```

OUTPUT
--

1. A SAM/BAM file that can be used for any downstream analyses.
2. A tab delimited UniProt report file.

```
#FORMAT

Count	UniProtKB	ID	Organism	Protein Names	Genes	Pathway	Features	Gene Ontology	Reviewd	Existence	Comments
```

- Count: the number of reads mapping to that UniProt entry
- UniProtKB: The ID containing the Gene short-code and species of origin
- ID: The Uniprot code
- Organims: The Organims from which the Uniprot ID is derived. Note that one should use this to generate a taxonomic profile of your sample
- Protein Names
- Genes
- Pathway	Features
- Gene Ontology
- Reviewd
- Existence
- Comments


[![PALADIN Wiki](https://github.com/twestbrookunh/paladin/wiki)]

[![Join the chat at https://gitter.im/twestbrookunh/paladin](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/twestbrookunh/paladin?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


