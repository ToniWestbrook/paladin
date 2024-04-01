#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include "align.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "bwa.h"
#include "protein.h"
#include "translations.h"
#include "uniprot.h"

static void *process(void *shared, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	const mem_opt_t *opt = aux->opt;
	int i;

	if (step == 0) {
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		if (!aux->copy_comment)
			for (i = 0; i < ret->n_seqs; ++i) {
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;

		logMessage(__func__, LOG_LEVEL_MESSAGE, "Read %d protein sequences (%ld AA)...\n", ret->n_seqs, (long)size);

		return ret;
	} else if (step == 1) {
		const bwaidx_t *idx = aux->idx;
		if (opt->flag & MEM_F_SMARTPE) {
			bseq1_t *sep[2];
			int n_sep[2];
			mem_opt_t tmp_opt = *opt;
			bseq_classify(data->n_seqs, data->seqs, n_sep, sep);

			logMessage(__func__, LOG_LEVEL_MESSAGE, "%d single-end sequences; %d paired-end sequences\n",  n_sep[0], n_sep[1]);

			if (n_sep[0]) {
				tmp_opt.flag &= ~MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, n_sep[0], sep[0], 0);
				for (i = 0; i < n_sep[0]; ++i)
					data->seqs[sep[0][i].id].sam = sep[0][i].sam;
			}
			if (n_sep[1]) {
				tmp_opt.flag |= MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed + n_sep[0], n_sep[1], sep[1], aux->pes0);
				for (i = 0; i < n_sep[1]; ++i)
					data->seqs[sep[1][i].id].sam = sep[1][i].sam;
			}
			free(sep[0]); free(sep[1]);
		} else mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, data->n_seqs, data->seqs, aux->pes0);
		aux->n_processed += data->n_seqs;

		return data;
	} else if (step == 2) {
		for (i = 0; i < data->n_seqs; ++i) {
			if (data->seqs[i].sam) err_fputs(data->seqs[i].sam, opt->outputStream);
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual); free(data->seqs[i].sam);
		}
		free(data->seqs); free(data);

		return 0;
	}
	return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0) {
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

int command_align(int argc, char *argv[]) {
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, ignore_alt = 0, no_mt_io = 0;
	int fixed_chunk_size = -1;
	gzFile fp, fp2 = 0;
	char *p, *rg_line = 0, *hdr_line = 0;
	const char *mode = 0;
	void *ko = 0, *ko2 = 0;
	mem_pestat_t pes[VALUE_DOMAIN];
	ktp_aux_t aux;
	FILE * reportPriStream = 0, * reportSecStream = 0;
	char * readsProName = 0, * indexProName = 0, * prefixName = 0;
	char * samName = 0, * reportPriName = 0, * reportSecName = 0;
    const char * proxyAddress;

	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, VALUE_DOMAIN * sizeof(mem_pestat_t));
	for (i = 0; i < VALUE_DOMAIN; ++i) pes[i].failed = 1;

	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
    proxyAddress = NULL;

	while ((c = getopt(argc, argv, "1epqabgnMCSVYJjf:F:z:u:k:o:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:P:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == 'u') opt->outputType = atoi(optarg);
		else if (c == 'f') opt->min_orf_len = atoi(optarg);
		else if (c == 'F') opt->min_orf_percent = atof(optarg);
        else if (c == 'z') opt->translations = convertTransArgs(optarg);
		else if (c == 'o') prefixName = optarg;
		else if (c == '1') no_mt_io = 1;
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		//else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		//else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'e') opt->flag |= MEM_F_SELF_OVLP;
		else if (c == 'F') opt->flag |= MEM_F_ALN_REG;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == 'b') opt->proteinFlag &= ~ALIGN_FLAG_BRUTE_ORF;
		else if (c == 'g') opt->proteinFlag |= ALIGN_FLAG_GEN_NT;
		else if (c == 'n') opt->proteinFlag |= ALIGN_FLAG_KEEP_PRO;
		else if (c == 'J') opt->proteinFlag &= ~ALIGN_FLAG_ADJUST_ORF;
        else if (c == 'p') opt->proteinFlag |= ALIGN_FLAG_MANUAL_PRO;
        else if (c == 'q') opt->proteinFlag |= ALIGN_FLAG_MANUAL_TRANSCRIPT;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') aux.copy_comment = 1;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
		else if (c == 'X') opt->mask_level = atof(optarg);
        else if (c == 'P') proxyAddress = optarg;
		else if (c == 'h') {
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		} else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		} else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 'H') {
			if (optarg[0] != '@') {
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0) {
					char *buf;
					buf = calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp)) {
						i = strlen(buf);
						assert(buf[i-1] == '\n'); // a long line
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		} else if (c == 'I') { // specify the insert size distribution
			aux.pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);

			logMessage(__func__, LOG_LEVEL_MESSAGE, "mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
			   									     pes[1].avg, pes[1].std, pes[1].high, pes[1].low);

		}
		else return 1;
	}

	if (rg_line) {
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 1 >= argc || optind + 3 < argc) {
		renderAlignUsage(opt);
		free(opt);
		return 1;
	}

	if (mode) {
		if (strcmp(mode, "intractg") == 0) {
			if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		} else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "pbread") == 0 || strcmp(mode, "ont2d") == 0) {
			if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "pbread") == 0) { // pacbio read-to-read setting; NOT working well!
				opt->flag |= MEM_F_ALL | MEM_F_SELF_OVLP | MEM_F_ALN_REG;
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.max_occ) opt->max_occ = 1000;
				if (!opt0.min_seed_len) opt->min_seed_len = 13;
				if (!opt0.max_chain_extend) opt->max_chain_extend = 25;
				if (opt0.drop_ratio == 0.) opt->drop_ratio = .001;
			} else if (strcmp(mode, "ont2d") == 0) {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			} else {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			logMessage(__func__, LOG_LEVEL_ERROR, "Unknown read type '%s'\n", mode);
			return 1; // FIXME memory leak
		}
	} else update_a(opt, &opt0);

    // Check for incompatible options
    if ((opt->proteinFlag & ALIGN_FLAG_MANUAL_PRO) && (opt->proteinFlag & ALIGN_FLAG_MANUAL_TRANSCRIPT)) {
		logMessage(__func__, LOG_LEVEL_ERROR, "Both protein and transcriptome mode selected.\n");
        return 1;
    }

    // Default to standard genetic code for translations
    if (!opt->translations) opt->translations = convertTransArgs("");

	// Create scoring weight matrix
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	// Prepare header check and filenames
	indexProName = malloc(strlen(argv[optind]) + 5);
	readsProName = malloc(strlen(argv[optind + 1]) + 5);
	sprintf(indexProName, "%s.pro", argv[optind]);
	sprintf(readsProName, "%s.pro", argv[optind + 1]);
	opt->indexInfo = getIndexHeader(indexProName);

	// Before loading index, ensure it's compatible
	switch (getIndexCompatible(opt->indexInfo)) {
		case INDEX_COMPATIBILITY_NONE:
			logMessage(__func__, LOG_LEVEL_ERROR,
					   "Index version %d.%d.%d is incompatible, please reindex the reference.\n",
						opt->indexInfo.version[0],
						opt->indexInfo.version[1],
						opt->indexInfo.version[2]);

			return 1;

		case INDEX_COMPATBILITY_FUTURE:
			logMessage(__func__, LOG_LEVEL_WARNING,
					   "Index version %d.%d.%d is newer than program version, compatibility not guaranteed.\n",
						opt->indexInfo.version[0],
						opt->indexInfo.version[1],
						opt->indexInfo.version[2]);

			break;
		default:
			break;
	}

	// Load index
	aux.idx = index_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		logMessage(__func__, LOG_LEVEL_MESSAGE, "Loading the index for reference '%s'...\n", argv[optind]);
		if ((aux.idx = index_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
	}
	else {
		logMessage(__func__, LOG_LEVEL_MESSAGE, "Loading the index from shared memory...\n");
	}

	if (ignore_alt)
		for (i = 0; i < aux.idx->bns->n_seqs; ++i)
			aux.idx->bns->anns[i].is_alt = 0;

	// Ready output files if requested (else stdout)
	if (prefixName != NULL) {
		if (opt->indexInfo.referenceType == 0) {
			logMessage(__func__, LOG_LEVEL_ERROR, "Reporting can only be used on prepared indices.\n");
 			return 1;
		}

		samName = malloc(strlen(prefixName) + 5);
		reportPriName = malloc(strlen(prefixName) + 23);
		reportSecName = malloc(strlen(prefixName) + 23);
		sprintf(samName, "%s.sam", prefixName);
		sprintf(reportPriName, "%s_uniprot.tsv", prefixName);

		if (opt->flag & MEM_F_ALL) {
			sprintf(reportPriName, "%s_uniprot_primary.tsv", prefixName);
			sprintf(reportSecName, "%s_uniprot_secondary.tsv", prefixName);
		}

		// Open files
		opt->outputStream = err_xopen_core(__func__, samName, "w");
		reportPriStream = err_xopen_core(__func__, reportPriName, "w");
		if (opt->flag & MEM_F_ALL) reportSecStream = err_xopen_core(__func__, reportSecName, "w");

	}

    if (opt->proteinFlag & ALIGN_FLAG_MANUAL_PRO) {
        // Protein input given, skip ORF detection
        sprintf(readsProName, "%s", argv[optind + 1]);
    }
    else {
        // Detect ORFs and write protein file
        writeReadsProtein(argv[optind + 1], readsProName, opt);
    }

	// Open ORFs sequence
	ko = kopen(readsProName, &fd);
	if (ko == 0) {
		logMessage(__func__, LOG_LEVEL_ERROR, "Failed to open file `%s'.\n", argv[optind + 1]);
		return 1;
	}

	fp = gzdopen(fd, "r");
	aux.ks = kseq_init(fp);
	if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			logMessage(__func__, LOG_LEVEL_WARNING, "When '-p' is in use, the second query file is ignored.\n");
		} else {
			ko2 = kopen(argv[optind + 2], &fd2);
			if (ko2 == 0) {
				logMessage(__func__, LOG_LEVEL_ERROR, "Failed to open file `%s'.\n", argv[optind + 2]);
				return 1;
			}
			fp2 = gzdopen(fd2, "r");
			aux.ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}

	// Render SAM header
	if (!(opt->flag & MEM_F_ALN_REG)) {
		bwa_print_sam_hdr(aux.idx->bns, hdr_line, opt->outputStream);
	}

	// Align and render
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
	kt_pipeline(no_mt_io? 1 : 2, process, &aux, 3);

	// Report number aligned
	renderNumberAligned(opt);

	// Generate UniProt report if requested
	if (prefixName != NULL) {
		renderUniprotReport(opt->outputType, 1, reportPriStream, proxyAddress);
		if (opt->flag & MEM_F_ALL) {
			renderUniprotReport(opt->outputType, 0, reportSecStream, proxyAddress);
		}
	}

	// Delete generated protein file unless requested otherwise
	if (!(opt->proteinFlag & ALIGN_FLAG_KEEP_PRO) && !(opt->proteinFlag & ALIGN_FLAG_MANUAL_PRO)) {
		remove(readsProName);
	}

	// Cleanup
	if (opt->outputStream != stdout) fclose (opt->outputStream);
	if (reportPriStream) fclose(reportPriStream);
	if (reportSecStream) fclose(reportSecStream);

	free(indexProName);
	free(readsProName);
	free(samName);
	free(reportPriName);
	free(reportSecName);
	free(hdr_line);
	free(opt);

	index_destroy(aux.idx);
	kseq_destroy(aux.ks);
	err_gzclose(fp); kclose(ko);
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		err_gzclose(fp2); kclose(ko2);
	}

	return 0;
}

int renderAlignUsage(const mem_opt_t * passOptions) {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: paladin align [options] <idxbase> <in.fq>\n\n");

	fprintf(stderr, "Gene detection options:\n\n");
    fprintf(stderr, "       -q            disable ORF detection and treat input as transcript sequence\n");
    fprintf(stderr, "       -p            disable ORF detection and treat input as protein sequence\n");
	fprintf(stderr, "       -b            disable brute force detection\n");
	fprintf(stderr, "       -J            do not adjust minimum ORF length (constant value) for shorter read lengths\n");
	fprintf(stderr, "       -f INT        minimum ORF length accepted (as constant value) [%d]\n", passOptions->min_orf_len);
	fprintf(stderr, "       -F FLOAT      minimum ORF length accepted (as percentage of read length) [%.2f]\n", passOptions->min_orf_percent);
    fprintf(stderr, "       -z INT[,...]  Genetic code used for translation (-z ? for full list) [1]\n");    

	fprintf(stderr, "\nAlignment options:\n\n");
	fprintf(stderr, "       -t INT        number of threads [%d]\n", passOptions->n_threads);
	fprintf(stderr, "       -k INT        minimum seed length [%d]\n", passOptions->min_seed_len);
	fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", passOptions->zdrop);
	fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", passOptions->split_factor);
	fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)passOptions->max_mem_intv);
    //fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", passOptions->split_width); // Disabled in BWA
	fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", passOptions->max_occ);
	fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", passOptions->drop_ratio);
	fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
	fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", passOptions->max_matesw);
	//fprintf(stderr, "       -S            skip mate rescue\n");
	//fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
	fprintf(stderr, "       -e            discard full-length exact matches\n");

	fprintf(stderr, "\nScoring options:\n\n");
	fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", passOptions->a);
	fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", passOptions->b);
	fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", passOptions->o_del, passOptions->o_ins);
	fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", passOptions->e_del, passOptions->e_ins);
	fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", passOptions->pen_clip5, passOptions->pen_clip3);
	fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", passOptions->pen_unpaired);
	fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
	fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
	fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
	fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
//		fprintf(stderr, "                     pbread: -k13 -W40 -c1000 -r10 -A1 -B1 -O1 -E1 -N25 -FeaD.001\n");
	fprintf(stderr, "\nInput/output options:\n\n");
	fprintf(stderr, "       -o STR        activate PALADIN reporting using STR as an output file prefix.  Files generated as follows:\n");
	fprintf(stderr, "                        STR.sam - alignment data (will not be sent to stdout)\n");
	fprintf(stderr, "                        STR_uniprot.tsv - Tab delimited UniProt report (normal alignment mode)\n");
	fprintf(stderr, "                        STR_uniprot_primary.tsv - Tab delimited UniProt report, primary alignments (all alignments mode)\n");
	fprintf(stderr, "                        STR_uniprot_secondary.tsv - Tab delimited UniProt report, secondary alignments (all alignments mode)\n\n");
	fprintf(stderr, "       -u INT        report type generated when using reporting and a UniProt reference [%d]\n", passOptions->outputType);
	fprintf(stderr, "                        0: Simple ID summary report\n");
	fprintf(stderr, "                        1: Detailed report (Contacts uniprot.org)\n\n");
    fprintf(stderr, "       -P STR        HTTP or SOCKS proxy address\n");
	fprintf(stderr, "       -g            generate detected ORF nucleotide sequence FASTA\n");
	fprintf(stderr, "       -n            keep protein sequence after alignment\n");
	//fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
	fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
	fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
	fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
	fprintf(stderr, "       -T INT        minimum score to output [%d]\n", passOptions->T);
	fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n", passOptions->max_XA_hits, passOptions->max_XA_hits_alt);
	fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
	fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
	fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
	fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
	fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
	fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
	fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
	fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
	fprintf(stderr, "                     FR orientation only. [inferred]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
	fprintf(stderr, "\n");

	return 1;
}
