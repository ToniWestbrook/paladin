#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include "utils.h"
#include "bwt.h"
#include "kvec.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

// Obtain the starting address of the interval in the BWT (post-interleave) for the specified index
uint32_t * getOccInterval(const bwt_t * passBWT, int64_t passSeqIdx) {
	int64_t numIntervals;

	numIntervals = passSeqIdx / 128;

	// Each interval stores 64-bit occurrences for each of the values + 128 packed values between
	return passBWT->bwt + (numIntervals * 2 * VALUE_DOMAIN) + (numIntervals * (128 / 4));
}

// Obtain unpacked byte value from the BWT (post-interleave) at the specified index
ubyte_t unpackBWTValue(const bwt_t * passBWT, int64_t passSeqIdx) {
	int64_t packShift;
	uint32_t * bwtValueStart;

	packShift = (8 * (sizeof(uint32_t) - 1)) - (8 * (passSeqIdx % sizeof(uint32_t)));

	bwtValueStart = getOccInterval(passBWT, passSeqIdx);
	bwtValueStart += sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN;

	return (bwtValueStart[(passSeqIdx % 128) / sizeof(uint32_t)] >> packShift) & 0xFF;
}

// Add to an occurrence array given a number of interleaved BWT occurrence byte values
void getOccPerWord(uint32_t passWord, bwtint_t retOcc[VALUE_DOMAIN], int64_t passCount) {
	ubyte_t * bytePtr;
	int64_t byteIdx;

	bytePtr = (ubyte_t *) & passWord;

	// Check each byte of word for occurrence of value
	for (byteIdx = (passCount - 1); byteIdx >= 0; byteIdx--) {
		retOcc[bytePtr[byteIdx]]++;
	}
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k) {
	bwtint_t x = k - (k > bwt->primary);

	x = unpackBWTValue(bwt, x);
	x = bwt->L2[x] + bwt_occ(bwt, k, x);

	return k == bwt->primary? 0 : x;
}

// bwt->bwt and bwt->occ must be precalculated
void bwt_cal_sa(bwt_t *bwt, int intv) {
	bwtint_t isa, sa, i; // S(isa) = sa
	int intv_round = intv;

	// Sanity checking
	kv_roundup32(intv_round);
	xassert(intv_round == intv, "SA sample interval is not a power of 2.");
	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	// Initialize Suffix Array
	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));

	// Calculate Suffix Array Value
	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % intv == 0) bwt->sa[isa/intv] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;

	// SA[0] is currently set to sequence length - set to maximum value
	bwt->sa[0] = (bwtint_t)-1;
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t sa = 0, mask = bwt->sa_intv - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv];
}

static inline int64_t __occ_aux(uint64_t y, int c) {
	unsigned char * bytePtr;
	int64_t count, byteIdx;

	bytePtr = (unsigned char *) &y;
	count = 0;

	for (byteIdx = 0 ; byteIdx < sizeof(y) ; byteIdx++) {
		if (bytePtr[byteIdx] == c) count++;
	}

	return count;
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c) {
	bwtint_t n;
	uint32_t *p, *end;

	//test
	bwtint_t cnt[VALUE_DOMAIN];
	memset(cnt, VALUE_DOMAIN, sizeof(bwtint_t));
	bwt_occ4(bwt, k, cnt);
	return cnt[c];

	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (bwtint_t)(-1)) return 0;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// Retrieve Occurrences at k/OCC_INTERVAL, then advance to first BWT cell
	n = ((bwtint_t*)(p = getOccInterval(bwt, k)))[c];
	p += sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN;

	// calculate Occ up to the last k/32
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

	// calculate Occ
	n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31; // corrected for the masked bits

	return n;
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
{
	bwtint_t _k, _l;
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		*ok = bwt_occ(bwt, k, c);
		*ol = bwt_occ(bwt, l, c);
	} else {
		bwtint_t m, n, i, j;
		uint32_t *p;
		if (k >= bwt->primary) --k;
		if (l >= bwt->primary) --l;
		//n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
		n = ((bwtint_t*)(p = getOccInterval(bwt, k)))[c];
		p += sizeof(bwtint_t);
		// calculate *ok
		j = k >> 5 << 5;
		for (i = k/OCC_INTERVAL*OCC_INTERVAL; i < j; i += 32, p += 2)
			n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
		m = n;
		n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
		if (c == 0) n -= ~k&31; // corrected for the masked bits
		*ok = n;
		// calculate *ol
		j = l >> 5 << 5;
		for (; i < j; i += 32, p += 2)
			m += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
		m += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
		if (c == 0) m -= ~l&31; // corrected for the masked bits
		*ol = m;
	}
}

void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[VALUE_DOMAIN])
{
	uint32_t *p, *end;

	// If K is set to max, set counts to 0
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, VALUE_DOMAIN * sizeof(bwtint_t));
		return;
	}

	k -= (k >= bwt->primary); // because $ is not in bwt

	// Copy current occurrence interval
	p = getOccInterval(bwt, k);
	memcpy(cnt, p, VALUE_DOMAIN * sizeof(bwtint_t));

	// Loop from end of interval to K in unpacked space
	p += sizeof(bwtint_t) / sizeof(uint32_t) * VALUE_DOMAIN;
	end = p + (((k / 16) - ((k&~OCC_INTV_MASK) / 16)) * 4) + ((k % 16) - (k % 4))/4;
	int supercount = 0;

	for (; p < end; ++p) {
		cnt[unpackBWTValue(bwt, k/128*128 + supercount)]++;
		cnt[unpackBWTValue(bwt, k/128*128 + supercount+1)]++;
		cnt[unpackBWTValue(bwt, k/128*128 + supercount+2)]++;
		cnt[unpackBWTValue(bwt, k/128*128 + supercount+3)]++;
		supercount+=4;
	}

	int testIdx;

	for (testIdx = 0 ; testIdx < (k % 4) + 1 ; testIdx++) {
		cnt[unpackBWTValue(bwt, k/128*128 + supercount++)]++;
	}
}

// STEP 2
// an analogy to bwt_occ4() but more efficient, requiring k <= l
void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[VALUE_DOMAIN], bwtint_t cntl[VALUE_DOMAIN])
{
	//bwtint_t _k, _l;
	//_k = k - (k >= bwt->primary);
	//_l = l - (l >= bwt->primary);

	bwt_occ4(bwt, k, cntk);
	bwt_occ4(bwt, l, cntl);

    // Temporarily bypassing the shortcut method for k <= l and calling the AA version of bwt_occ4 until (if possible) a faster method can be devised
/*
	if (_l>>OCC_INTV_SHIFT != _k>>OCC_INTV_SHIFT || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		bwtint_t x, y;
		uint32_t *p, tmp, *endk, *endl;
		k -= (k >= bwt->primary); // because $ is not in bwt
		l -= (l >= bwt->primary);
		p = getOccInterval(bwt, k);
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
		// prepare cntk[]
		endk = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4));
		endl = p + ((l>>4) - ((l&~OCC_INTV_MASK)>>4));
		for (x = 0; p < endk; ++p) x += __occ_aux4(bwt, *p);
		y = x;
		tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
		x += __occ_aux4(bwt, tmp) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		for (; p < endl; ++p) y += __occ_aux4(bwt, *p);
		tmp = *p & ~((1U<<((~l&15)<<1)) - 1);
		y += __occ_aux4(bwt, tmp) - (~l&15);
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
	*/
}

int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t k, l, ok, ol;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) break; // no match
	}
	if (k > l) return 0; // no match
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
	int i;
	bwtint_t k, l, ok, ol;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}

/*********************
 * Bidirectional BWT *
 *********************/
// STEP 1
void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[VALUE_DOMAIN], int is_back)
{
	bwtint_t tk[VALUE_DOMAIN], tl[VALUE_DOMAIN];
	int i, copyIdx;

	// Accumulate occurrences
	bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != VALUE_DOMAIN; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}

	// Copy values returned from occurrence accumulators
	ok[VALUE_DOMAIN - 1].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);

	for (copyIdx = VALUE_DOMAIN - 2 ; copyIdx >= 0 ; copyIdx--) {
		ok[copyIdx].x[is_back] = ok[copyIdx + 1].x[is_back] + ok[copyIdx + 1].x[2];
	}
}

static void bwt_reverse_intvs(bwtintv_v *p)
{
	if (p->n > 1) {
		int j;
		for (j = 0; j < p->n>>1; ++j) {
			bwtintv_t tmp = p->a[p->n - 1 - j];
			p->a[p->n - 1 - j] = p->a[j];
			p->a[j] = tmp;
		}
	}
}
// NOTE: $max_intv is not currently used in BWA-MEM
int bwt_smem1a(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	int i, j, c, ret;
	bwtintv_t ik, ok[VALUE_DOMAIN];
	bwtintv_v a[2], *prev, *curr, *swap;

	mem->n = 0;
	if (q[x] >= VALUE_DEFINED) return x + 1;
	//if (q[x] > 3) return x + 1;
	if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
	kv_init(a[0]); kv_init(a[1]);
	prev = tmpvec && tmpvec[0]? tmpvec[0] : &a[0]; // use the temporary vector if provided
	curr = tmpvec && tmpvec[1]? tmpvec[1] : &a[1];

	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base

	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (ik.x[2] < max_intv) { // an interval small enough
			kv_push(bwtintv_t, *curr, ik);
			break;
		} else if (q[i] < VALUE_DEFINED) { // an amino acid
			c = VALUE_DEFINED - 1 - q[i]; // complement of q[i]
			bwt_extend(bwt, &ik, ok, 0);

			if (ok[c].x[2] != ik.x[2]) { // change of the interval size
				kv_push(bwtintv_t, *curr, ik);
				if (ok[c].x[2] < min_intv) break; // the interval size is too small to be extended further
			}
			ik = ok[c]; ik.info = i + 1;
		} else { // an ambiguous base
			kv_push(bwtintv_t, *curr, ik);
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}

	if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
	bwt_reverse_intvs(curr); // s.t. smaller intervals (i.e. longer matches) visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? -1 : q[i] < VALUE_DEFINED? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
		//c = i < 0? -1 : q[i] < 4? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			bwtintv_t *p = &prev->a[j];
			if (c >= 0 && ik.x[2] >= max_intv) bwt_extend(bwt, p, ok, 1);
			if (c < 0 || ik.x[2] < max_intv || ok[c].x[2] < min_intv) { // keep the hit if reaching the beginning or an ambiguous base or the intv is small enough
				if (curr->n == 0) { // test curr->n>0 to make sure there are no longer matches
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
						kv_push(bwtintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			} else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]) {
				ok[c].info = p->info;
				kv_push(bwtintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}

	bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

	if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
	if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
	return ret;
}

int bwt_smem1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	return bwt_smem1a(bwt, len, q, x, min_intv, 0, mem, tmpvec);
}

int bwt_seed_strategy1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem)
{
	int i, c;
	bwtintv_t ik, ok[VALUE_DOMAIN];

	memset(mem, 0, sizeof(bwtintv_t));
	if (q[x] >= VALUE_DEFINED) return x + 1;
	//if (q[x] > 3) return x + 1;
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	for (i = x + 1; i < len; ++i) { // forward search
		if (q[i] < VALUE_DEFINED) { // an amino acid
		c = VALUE_DEFINED - 1 - q[i]; // complement of q[i]
		bwt_extend(bwt, &ik, ok, 0);
		if (ok[c].x[2] < max_intv && i - x >= min_len) {
			*mem = ok[c];
			mem->info = (uint64_t)x<<32 | (i + 1);
			return i + 1;
		}
		ik = ok[c];
		} else return i + 1;
	}
	return len;
}

/*************************
 * Read/write BWT and SA *
 *************************/

void bwt_dump_bwt(const char *fn, const bwt_t *bwt) {
	FILE *fp;

	fp = xopen(fn, "wb");

	// Format: <primary><L2[1::]><bwt>
	// Sizes: 8, Domain * 8, BWTSize * 4
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), VALUE_DOMAIN, fp);
	err_fwrite(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp);
	err_fflush(fp);
	err_fclose(fp);
}

void bwt_dump_sa(const char *fn, const bwt_t *bwt) {
	FILE *fp;

	fp = xopen(fn, "wb");

	// Format: <primary><L2[1::]><SA Interval><Seq Length><SA[1::]>
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), VALUE_DOMAIN, fp);
	err_fwrite(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fwrite(&bwt->seq_len, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	err_fflush(fp);
	err_fclose(fp);
}

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

void bwt_restore_sa(const char *fn, bwt_t *bwt) {
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = xopen(fn, "rb");

	// Format: <primary><L2[1::]><SA Interval><Seq Length><SA[1::]>
	err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	err_fread_noeof(skipped, sizeof(bwtint_t), VALUE_DOMAIN, fp); // skip
	err_fread_noeof(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	fread_fix(fp, sizeof(bwtint_t) * (bwt->n_sa - 1), bwt->sa + 1);
	err_fclose(fp);
}

bwt_t *bwt_restore_bwt(const char *fn) {
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = xopen(fn, "rb");
	err_fseek(fp, 0, SEEK_END);

	// Populate BWT Size and allocate BWT
	bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * (VALUE_DOMAIN + 1)) / sizeof(uint32_t);
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, sizeof(uint32_t));

	// Populate Primary, L2[1::], and BWT
	err_fseek(fp, 0, SEEK_SET);
	err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), VALUE_DOMAIN, fp);
	fread_fix(fp, bwt->bwt_size * sizeof(uint32_t), bwt->bwt);

	// Sequence length in end of L2 array
	bwt->seq_len = bwt->L2[VALUE_DOMAIN];

	err_fclose(fp);

	return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}
