#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);
}



/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param w_min  minimum strobe offset
 * @param w_max  maximum strobe offset
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch_altstrobes(void *km, const char *str, int len, int w, int k, int k_min, int w_min, int w_max, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, buf_temp, buf_temp2, min_pos, kmer_span, strobe_pos_next, w_buf_pos, w_poss, counter = 0;
	mm128_t buf_fw[256], min = { UINT64_MAX, UINT64_MAX };
	mm128_t buf_rev[256];
	tiny_queue_t tq;
	int w_min1 = w_min - k/2;
	int w_min2 = w_min + k/2;
	int w_max1 = w_max - k/2;
	int w_max2 = w_max + k/2;
	int w_max_w = 2*w_max2 + w;
	uint64_t min_val = UINT64_MAX;
	uint64_t info_2k = UINT64_MAX;
	int strobe_span = 0;
	int strobe_len = w_max2+k-1;
	int min_start = k+w_max2+w;

	assert(len > 0 && (w_max_w > 0 && w_max_w < 256) && (k > 0 && k <= 9)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf_fw, 0xff, w_max_w * 16);
	memset(buf_rev, 0xff, w_max_w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info_fw = { UINT64_MAX, UINT64_MAX };
		mm128_t info_rev = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			// if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			// z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info_fw.x = hash64(kmer[0], mask); // << 8 | kmer_span;
				// info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
				info_rev.x = hash64(kmer[1], mask); // << 8 | kmer_span;
				// info_rev.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf_fw[buf_pos] = info_fw; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		buf_rev[buf_pos] = info_rev; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		// fprintf(stderr, "L: %d, Buf_pos: %d\n", l, buf_pos);
	
		if (l >= min_start){
			// fprintf(stderr, "Buf_pos: %d, w_buf %d, Min: %d, Fwd: %d, Rev: %d\n", buf_pos, w_buf_pos, min.x, buf_fw[(buf_pos+1+w)%w_max_w].x, info_rev.x);
			buf_temp = (w_max_w+buf_pos-w_max2)%w_max_w;
			if ((buf_fw[buf_temp].x <= min.x) && (buf_fw[buf_temp].x <= buf_rev[buf_temp].x)){  // a new minimum
				min_val = UINT64_MAX;
				if (buf_fw[buf_temp].x % 2 == 0){ // sample k-2k
					for (j = w_min1+buf_temp; j < w_max1+buf_temp; j++){
						uint64_t res = buf_fw[buf_temp].x ^ buf_fw[j%w_max_w].x ^ buf_fw[(j+k)%w_max_w].x;
						if (res < min_val){
							min_val = res;
							strobe_pos_next = j;
						}
					}
					// fprintf(stderr, "[Fwd k-2k] w_buf: 0, %d, %d, %d, %d, %d, %d\n", (info_fw.x << (2*k)) ^ (buf_fw[strobe_pos_next%w_max_w].x << k) ^ (buf_fw[(strobe_pos_next+k)%w_max_w].x),info_fw.x,buf_fw[strobe_pos_next%w_max_w].x,buf_fw[(strobe_pos_next+k)%w_max_w].x,(i-w_max2+strobe_span-k), strobe_span);
					info_fw.x = ((buf_fw[buf_temp].x << (2*k)) ^ (buf_fw[strobe_pos_next%w_max_w].x << k) ^ (buf_fw[(strobe_pos_next+k)%w_max_w].x)) << 8 | strobe_len;
					info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
					// fprintf(stderr, "Fwd_k-2k Normal: %d \n", i);
					kv_push(mm128_t, km, *p, info_fw);
				} else { // sample 2k-k
					info_2k = buf_fw[buf_temp].x ^ buf_fw[(buf_temp+k)%w_max_w].x;
					for (j = w_min2+buf_temp; j < w_max2+buf_temp; j++){
						uint64_t res = (info_2k ^ buf_fw[j%w_max_w].x);
						if (res < min_val){
							min_val = res;
							strobe_pos_next = j;
						}
					}
					// fprintf(stderr, "[Fwd 2k-k] w_buf: 0, %d, %d, %d, %d, %d, %d\n", (info_fw.x << (2*k)) ^ (buf_fw[(buf_temp+k)%w_max_w].x << k) ^ (buf_fw[strobe_pos_next%w_max_w].x), info_fw.x, buf_fw[(buf_temp+k)%w_max_w].x, buf_fw[strobe_pos_next%w_max_w].x, (i-w_max2+strobe_span-k), strobe_span);
					info_fw.x = ((buf_fw[buf_temp].x << (2*k)) ^ (buf_fw[(buf_temp+k)%w_max_w].x << k) ^ (buf_fw[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
					info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
					// fprintf(stderr, "Fwd_2k-k Normal: %d \n", i);
					kv_push(mm128_t, km, *p, info_fw);
				}
				w_buf_pos = 0;
				min = buf_fw[buf_temp];
			}
			else if ((buf_rev[buf_temp].x <= min.x) && (l >= min_start+w_max2)) {
				min_val = UINT64_MAX;
				if (buf_rev[buf_temp].x % 2 == 0){ // sample k-2k
					for (j = w_max_w+buf_temp-w_min1; j > w_max_w+buf_temp-w_max1; j--){
						uint64_t res = (buf_rev[buf_temp].x ^ buf_rev[j%w_max_w].x ^ buf_rev[(j-k)%w_max_w].x);
						if (res < min_val){
							min_val = res;
							strobe_pos_next = j;
						}
					}
					// fprintf(stderr, "[Rev k-2k] w_buf: 0, %d, %d, %d, %d, %d, %d\n", (info_rev.x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x << k) ^ (buf_rev[(strobe_pos_next-k)%w_max_w].x), info_rev.x, buf_rev[strobe_pos_next%w_max_w].x, buf_rev[(strobe_pos_next-k)%w_max_w].x, i, (-strobe_pos_next+w_max_w+buf_pos+2*k));
					info_rev.x = ((buf_rev[buf_temp].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x << k) ^ (buf_rev[(strobe_pos_next-k)%w_max_w].x))  << 8 | strobe_len;
					info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max2)<<1 | 1;
					// fprintf(stderr, "Rev_k-2k Normal: %d \n", i);
					kv_push(mm128_t, km, *p, info_rev);
				} else { // sample 2k-k
					info_2k = buf_rev[buf_temp].x ^ buf_rev[(w_max_w+buf_temp-k)%w_max_w].x;
					for (j = w_max_w+buf_temp-w_min2; j > w_max_w+buf_temp-w_max2; j--){
						uint64_t res = (info_2k ^ buf_rev[j%w_max_w].x);
						if (res < min_val){
							min_val = res;
							strobe_pos_next = j;
						}
					}
					// fprintf(stderr, "[Rev 2k-k] w_buf: 0, %d, %d, %d, %d, %d, %d\n", (info_rev.x << (2*k)) ^ (buf_rev[(w_max_w+buf_pos-k)%w_max_w].x << k) ^ (buf_rev[strobe_pos_next%w_max_w].x), info_rev.x, buf_rev[(w_max_w+buf_pos-k)%w_max_w].x, buf_rev[strobe_pos_next%w_max_w].x, i, (-strobe_pos_next+w_max_w+buf_pos+k));
					info_rev.x = ((buf_rev[buf_temp].x << (2*k)) ^ (buf_rev[(w_max_w+buf_temp-k)%w_max_w].x << k) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
					info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max2)<<1 | 1;
					// fprintf(stderr, "Rev_2k-k Normal: %d \n", i);
					kv_push(mm128_t, km, *p, info_rev);
				}
				w_buf_pos = 0;
				min = buf_rev[buf_temp];
			}
			if (++w_buf_pos == w){  // old min has moved outside the window
				// fprintf(stderr, "Moved outside window\n");
				min.x =  UINT64_MAX;
				buf_temp = w_max_w+buf_pos-w_max2;
				for (w_poss = (w-2); w_poss >= 0; --w_poss){
					// fprintf(stderr, "\tw_poss: %d, Fwd: %d, Rev: %d\n", w_poss, buf_fw[(buf_temp-w_poss)%w_max_w].x, buf_rev[(w_max_w+buf_pos-w_poss)%w_max_w].x);
					if (buf_fw[(buf_temp-w_poss)%w_max_w].x < min.x){
						w_buf_pos = w_poss;
						min = buf_fw[(buf_temp-w_poss)%w_max_w];
					}
					if (buf_rev[(buf_temp-w_poss)%w_max_w].x < min.x){
						w_buf_pos = w_poss;
						min = buf_rev[(buf_temp-w_poss)%w_max_w];
					}
				}
				for (w_poss = (w-2); w_poss >= 0; --w_poss){
					// fprintf(stderr, "w: %d, w_poss: %d\n", w, w_poss);
					if (buf_fw[(buf_temp-w_poss)%w_max_w].x == min.x){
						buf_temp2 = (buf_temp-w_poss)%w_max_w;
						min_val = UINT64_MAX;
						if (buf_fw[buf_temp2].x % 2 == 0){ // sample k-2k
							for (j = w_min1+buf_temp2; j < w_max1+buf_temp2; j++){
								uint64_t res = buf_fw[buf_temp2].x ^ buf_fw[j%w_max_w].x ^ buf_fw[(j+k)%w_max_w].x;
								if (res < min_val){
									min_val = res;
									strobe_pos_next = j;
								}
							}
							// fprintf(stderr, "[Fwd k-2k] w_buf: %d, %d, %d, %d, %d, %d, %d\n", w_buf_pos, (info_fw.x << (2*k)) ^ (buf_fw[strobe_pos_next%w_max_w].x << k) ^ (buf_fw[(strobe_pos_next+k)%w_max_w].x),info_fw.x,buf_fw[strobe_pos_next%w_max_w].x,buf_fw[(strobe_pos_next+k)%w_max_w].x,(i-w_max_w+w+strobe_span-k), strobe_span);
							info_fw.x = ((buf_fw[buf_temp2].x << (2*k)) ^ (buf_fw[strobe_pos_next%w_max_w].x << k) ^ (buf_fw[(strobe_pos_next+k)%w_max_w].x)) << 8 | strobe_len;
							info_fw.y = (uint64_t)rid<<32 | (uint32_t)(i-w_poss)<<1 | 0;
							// fprintf(stderr, "Fwd_k-2k Window: %d %d \n", i, w_poss);
							kv_push(mm128_t, km, *p, info_fw);
						} else { // sample 2k-k
							info_2k = buf_fw[(buf_temp2)%w_max_w].x ^ buf_fw[(buf_temp2+k)%w_max_w].x;
							for (j = w_min2+buf_temp2; j < w_max2+buf_temp2; j++){
								uint64_t res = (info_2k ^ buf_fw[j%w_max_w].x);
								if (res < min_val){
									min_val = res;
									strobe_pos_next = j;
								}
							}
							// fprintf(stderr, "[Fwd 2k-k] w_buf: %d, %d, %d, %d, %d, %d, %d\n", w_buf_pos, (info_fw.x << (2*k)) ^ (buf_fw[(buf_temp2+k)%w_max_w].x << k) ^ (buf_fw[strobe_pos_next%w_max_w].x), info_fw.x, buf_fw[(buf_temp2+k)%w_max_w].x, buf_fw[strobe_pos_next%w_max_w].x, (i-w_max_w+w+strobe_span-k), strobe_span);
							info_fw.x = ((buf_fw[buf_temp2].x << (2*k)) ^ (buf_fw[(buf_temp2+k)%w_max_w].x << k) ^ (buf_fw[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
							info_fw.y = (uint64_t)rid<<32 | (uint32_t)(i-w_poss)<<1 | 0;
							// fprintf(stderr, "Fwd_2k-k Window: %d %d\n", i, w_poss);
							kv_push(mm128_t, km, *p, info_fw);
						}
					}

					if ((buf_rev[(buf_temp-w_poss)%w_max_w].x == min.x) && (l >= min_start+w_max2)){
						buf_temp2 = (buf_temp-w_poss)%w_max_w;
						min_val = UINT64_MAX;
						if (buf_rev[buf_temp2].x % 2 == 0){ // sample k-2k
							for (j = w_max_w+buf_temp2-w_min1; j > w_max_w+buf_temp2-w_max1; j--){
								uint64_t res = (buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x ^ buf_rev[(j-k)%w_max_w].x);
								if (res < min_val){
									min_val = res;
									strobe_pos_next = j;
								}
							}
							// fprintf(stderr, "[Rev k-2k] w_buf: %d, %d, %d, %d, %d, %d, %d\n", w_buf_pos, (buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x << k) ^ (buf_rev[(strobe_pos_next-k)%w_max_w].x), buf_rev[buf_temp2].x, buf_rev[strobe_pos_next%w_max_w].x, buf_rev[(strobe_pos_next-k)%w_max_w].x, (i+w-w_poss), (-strobe_pos_next+w_max_w+buf_temp2+2*k));
							info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x << k) ^ (buf_rev[(strobe_pos_next-k)%w_max_w].x))  << 8 | strobe_len;
							info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max2-w_poss)<<1 | 1;
							// fprintf(stderr, "Rev_k-2k Window: %d %d \n", i, w_poss);
							kv_push(mm128_t, km, *p, info_rev);
						} else { // sample 2k-k
							info_2k = buf_rev[buf_temp2].x ^ buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x;
							for (j = w_max_w+buf_temp2-w_min2; j > w_max_w+buf_temp2-w_max2; j--){
								uint64_t res = (info_2k ^ buf_rev[j%w_max_w].x);
								if (res < min_val){
									min_val = res;
									strobe_pos_next = j;
								}
							}
							// fprintf(stderr, "[Rev 2k-k] w_buf: %d, %d, %d, %d, %d, %d, %d\n", w_buf_pos, (buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[(w_max_w+buf_pos-k)%w_max_w].x << k) ^ (buf_rev[strobe_pos_next%w_max_w].x), buf_rev[buf_temp2].x, buf_rev[(w_max_w+buf_pos-k)%w_max_w].x, buf_rev[strobe_pos_next%w_max_w].x, (i+w-w_poss), (-strobe_pos_next+w_max_w+buf_temp2+k));
							info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x << k) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
							info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max2-w_poss)<<1 | 1;
							// fprintf(stderr, "Rev_2k-k Window: %d %d \n", i, w_poss);
							kv_push(mm128_t, km, *p, info_rev);
						}
					}
				}
				++w_buf_pos;
			}
		}
		if (++buf_pos == w_max_w) buf_pos = 0;
	}
  // compute rev minimizers on last w_max pos:
	// fprintf(stderr, "compute rev minimzer\n");
	mm128_t info_rev = { UINT64_MAX, UINT64_MAX };
	counter = 0;
	for (buf_temp = (w_max_w+buf_pos-w_max); buf_temp < (w_max_w+buf_pos); ++buf_temp){
		++counter;
		buf_temp2 = buf_temp%w_max_w;
		if (buf_rev[buf_temp2].x <= min.x) {
			min_val = UINT64_MAX;
			if (buf_rev[buf_temp2].x % 2 == 0){ // sample k-2k
				for (j = w_max_w+buf_temp2-w_min1; j > w_max_w+buf_temp2-w_max1; j--){
					uint64_t res = (buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x ^ buf_rev[(j-k)%w_max_w].x);
					if (res < min_val){
						min_val = res;
						strobe_pos_next = j;
					}
				}
				// fprintf(stderr, "[Rev k-2k] w_buf: 0, %d, %d, %d, %d, %d, %d\n", (info_rev.x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x << k) ^ (buf_rev[(strobe_pos_next-k)%w_max_w].x), info_rev.x, buf_rev[strobe_pos_next%w_max_w].x, buf_rev[(strobe_pos_next-k)%w_max_w].x, i, (-strobe_pos_next+w_max_w+buf_pos+2*k));
				info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x << k) ^ (buf_rev[(strobe_pos_next-k)%w_max_w].x))  << 8 | strobe_len;
				info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max2+counter-1)<<1 | 1;
				// fprintf(stderr, "Rev_k-2k End: %d %d\n", i, counter);
				kv_push(mm128_t, km, *p, info_rev);
			} else { // sample 2k-k
				info_2k = buf_rev[buf_temp2].x ^ buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x;
				for (j = w_max_w+buf_temp2-w_min2; j > w_max_w+buf_temp2-w_max2; j--){
					uint64_t res = (info_2k ^ buf_rev[j%w_max_w].x);
					if (res < min_val){
						min_val = res;
						strobe_pos_next = j;
					}
				}
				// fprintf(stderr, "[Rev 2k-k] w_buf: 0, %d, %d, %d, %d, %d, %d\n", (info_rev.x << (2*k)) ^ (buf_rev[(w_max_w+buf_pos-k)%w_max_w].x << k) ^ (buf_rev[strobe_pos_next%w_max_w].x), info_rev.x, buf_rev[(w_max_w+buf_pos-k)%w_max_w].x, buf_rev[strobe_pos_next%w_max_w].x, i, (-strobe_pos_next+w_max_w+buf_pos+k));
				info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x << k) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
				info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max2+counter-1)<<1 | 1;
				// fprintf(stderr, "Rev_2k-k End: %d %d\n", i, counter);
				kv_push(mm128_t, km, *p, info_rev);
			}
		}
		if (++w_buf_pos == w){  // old min has moved outside the window
			// fprintf(stderr, "Moved outside window\n");
			min.x = UINT64_MAX;
			for (w_poss = (w-2); w_poss >= 0; --w_poss){
				// fprintf(stderr, "\tw_poss: %d, Fwd: %d, Rev: %d\n", w_poss, buf_fw[(buf_temp-w_poss)%w_max_w].x, buf_rev[(w_max_w+buf_pos-w_poss)%w_max_w].x);
				if (buf_rev[(buf_temp-w_poss)%w_max_w].x < min.x){
					w_buf_pos = w_poss;
					min = buf_rev[(buf_temp-w_poss)%w_max_w];
				}
			}
			for (w_poss = (w-2); w_poss >= 0; --w_poss){
				if (buf_rev[(buf_temp-w_poss)%w_max_w].x == min.x) {
					min_val = UINT64_MAX;
					buf_temp2 = (buf_temp - w_poss)%w_max_w;
					if (buf_rev[buf_temp2].x % 2 == 0){ // sample k-2k
						for (j = w_max_w+buf_temp2-w_min1; j > w_max_w+buf_temp2-w_max1; j--){
							uint64_t res = (buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x ^ buf_rev[(j-k)%w_max_w].x);
							if (res < min_val){
								min_val = res;
								strobe_pos_next = j;
							}
						}
						// fprintf(stderr, "[Rev k-2k] w_buf: %d, %d, %d, %d, %d, %d, %d\n", w_buf_pos, (buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x << k) ^ (buf_rev[(strobe_pos_next-k)%w_max_w].x), buf_rev[buf_temp2].x, buf_rev[strobe_pos_next%w_max_w].x, buf_rev[(strobe_pos_next-k)%w_max_w].x, (i+w-w_poss), (-strobe_pos_next+w_max_w+buf_temp2+2*k));
						info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x << k) ^ (buf_rev[(strobe_pos_next-k)%w_max_w].x))  << 8 | strobe_len;
						info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max2-w_poss+counter-1)<<1 | 1;
						// fprintf(stderr, "Rev_k-2k EndWindow: %d %d %d\n", i, w_poss, counter);
						kv_push(mm128_t, km, *p, info_rev);
					} else { // sample 2k-k
						info_2k = buf_rev[buf_temp2].x ^ buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x;
						for (j = w_max_w+buf_temp2-w_min2; j > w_max_w+buf_temp2-w_max2; j--){
							uint64_t res = (info_2k ^ buf_rev[j%w_max_w].x);
							if (res < min_val){
								min_val = res;
								strobe_pos_next = j;
							}
						}
						// fprintf(stderr, "[Rev 2k-k] w_buf: %d, %d, %d, %d, %d, %d, %d\n", w_buf_pos, (buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[(w_max_w+buf_pos-k)%w_max_w].x << k) ^ (buf_rev[strobe_pos_next%w_max_w].x), buf_rev[buf_temp2].x, buf_rev[(w_max_w+buf_pos-k)%w_max_w].x, buf_rev[strobe_pos_next%w_max_w].x, (i+w-w_poss), (-strobe_pos_next+w_max_w+buf_temp2+k));
						info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x << k) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
						info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max2-w_poss+counter-1)<<1 | 1;
						// fprintf(stderr, "Rev_2k-k EndWindow: %d %d %d\n", i, w_poss, counter);
						kv_push(mm128_t, km, *p, info_rev);
					}
				}
			}
			++w_buf_pos;
		}
	}

}


/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param w_min  minimum strobe offset
 * @param w_max  maximum strobe offset
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch_mixedstrobes(void *km, const char *str, int len, int w, int k, int k_min, int w_min, int w_max, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, buf_temp, buf_temp2, min_pos, kmer_span, strobe_pos_next, w_buf_pos, w_poss, counter = 0;
	mm128_t buf_fw[256], min = { UINT64_MAX, UINT64_MAX };
	mm128_t buf_rev[256];
	tiny_queue_t tq;
	uint64_t min_val = UINT64_MAX;
	int w_max_w = 2*w_max + w;
	int strobe_span = 0;
	int strobe_len = w_max+k-1;
    int min_start = k+w_max+w;

	assert(len > 0 && (w_max_w > 0 && w_max_w < 256) && (k > 0 && k <= 14)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf_fw, 0xff, w_max_w * 16);
	memset(buf_rev, 0xff, w_max_w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info_fw = { UINT64_MAX, UINT64_MAX };
		mm128_t info_rev = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			// if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			// z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info_fw.x = hash64(kmer[0], mask); // << 8 | kmer_span;
				// info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
				info_rev.x = hash64(kmer[1], mask); // << 8 | kmer_span;
				// info_rev.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf_fw[buf_pos] = info_fw; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		buf_rev[buf_pos] = info_rev; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below

		if (l >= min_start){
			// fprintf(stderr, "Buf_pos: %d, w_buf %d, Min: %d, Fwd: %d, Rev: %d\n", buf_pos, w_buf_pos, min.x, buf_fw[(buf_pos+1+w)%w_max_w].x, info_rev.x);
			buf_temp = (w_max_w+buf_pos-w_max)%w_max_w;
			if ((buf_fw[buf_temp].x <= min.x) && (buf_fw[buf_temp].x <= buf_rev[buf_temp].x)){
				min_val = UINT64_MAX;
				if (buf_fw[buf_temp].x % 100 < 80){ // sample randstrobes
					info_fw.x = buf_fw[buf_temp].x;
					for (j = w_min+buf_temp; j < w_max+buf_temp; j++){
						uint64_t res = info_fw.x ^ buf_fw[j%w_max_w].x;
						if (res < min_val){
							min_val = res;
							strobe_pos_next = j;
						}
					}
					info_fw.x = ((info_fw.x << (2*k)) ^ (buf_fw[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
					info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
					kv_push(mm128_t, km, *p, info_fw);
				} else { // sample k-mers
					// fprintf(stderr, "Fwd: %d, %d\n", buf_fw[(buf_pos+1)%w_max_w].x, buf_fw[(buf_pos+1+k)%w_max_w].x);
					info_fw.x = ((buf_fw[buf_temp].x << (2*k)) ^ (buf_fw[(buf_temp+k)%w_max_w].x)) << 8 | strobe_len;
					info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
					kv_push(mm128_t, km, *p, info_fw);
				}
				w_buf_pos = 0;
				min = buf_fw[buf_temp];
			}
			else if ((buf_rev[buf_temp].x <= min.x) && (l >= min_start+w_max)) {
				min_val = UINT64_MAX;
				if (buf_rev[buf_temp].x % 100 < 80){ // sample randstrobes
					for (j = w_max_w+buf_temp-w_min; j > w_max_w+buf_temp-w_max; j--){
						uint64_t res = buf_rev[buf_temp].x ^ buf_rev[j%w_max_w].x;
						if (res < min_val){
							min_val = res;
							strobe_pos_next = j;
						}
					}
					info_rev.x = ((buf_rev[buf_temp].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
					info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max)<<1 | 1;
					kv_push(mm128_t, km, *p, info_rev);
				} else { // sample k-mers
					// fprintf(stderr, "Rev: %d, %d\n", buf_rev[(w_max_w+buf_temp-k)%w_max_w].x, buf_rev[buf_temp].x);
					info_rev.x = ((buf_rev[buf_temp].x << (2*k)) ^ (buf_rev[(w_max_w+buf_temp-k)%w_max_w].x)) << 8 | strobe_len;
					info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max)<<1 | 1;
					kv_push(mm128_t, km, *p, info_rev);
				}
				w_buf_pos = 0;
				min = buf_rev[buf_temp];
			}

			if (++w_buf_pos == w){  // old min has moved outside the window
				// fprintf(stderr, "Moved outside window\n");
				min.x = UINT64_MAX;
                buf_temp = w_max_w+buf_pos-w_max;
				for (w_poss = (w-2); w_poss >= 0; --w_poss){
					// fprintf(stderr, "\tw_poss: %d, Fwd: %d, Rev: %d\n", w_poss, buf_fw[(buf_temp-w_poss)%w_max_w].x, buf_rev[(w_max_w+buf_pos-w_poss)%w_max_w].x);
					if (buf_fw[(buf_temp-w_poss)%w_max_w].x < min.x){
						w_buf_pos = w_poss;
						min = buf_fw[(buf_temp-w_poss)%w_max_w];
					}
					if ((buf_rev[(buf_temp-w_poss)%w_max_w].x < min.x) && (l >= min_start+w_max)){
						w_buf_pos = w_poss;
						min = buf_rev[(buf_temp-w_poss)%w_max_w];
					}
				}
				for (w_poss = (w-2); w_poss >= 0; --w_poss){
					if (buf_fw[(buf_temp-w_poss)%w_max_w].x == min.x){
						buf_temp2 = (buf_temp - w_poss)%w_max_w;
						min_val = UINT64_MAX;
						if (buf_fw[buf_temp2].x % 100 < 80){ // sample randstrobes
							info_fw.x = buf_fw[buf_temp2].x;
							for (j = w_min+buf_temp2; j < w_max+buf_temp2; j++){
								uint64_t res = info_fw.x ^ buf_fw[j%w_max_w].x;
								if (res < min_val){
									min_val = res;
									strobe_pos_next = j;
								}
							}
							info_fw.x = ((info_fw.x << (2*k)) ^ (buf_fw[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
							info_fw.y = (uint64_t)rid<<32 | (uint32_t)(i-w_poss)<<1 | 0;
							kv_push(mm128_t, km, *p, info_fw);
						} else { // sample k-mers
							info_fw.x = ((buf_fw[buf_temp2].x << (2*k)) ^ (buf_fw[(buf_temp2+k)%w_max_w].x)) << 8 | strobe_len;
							info_fw.y = (uint64_t)rid<<32 | (uint32_t)(i-w_poss)<<1 | 0;
							kv_push(mm128_t, km, *p, info_fw);
						}
					}
					if ((buf_rev[(buf_temp-w_poss)%w_max_w].x == min.x) && (l >= min_start+w_max)) {
						min_val = UINT64_MAX;
						buf_temp2 = (buf_temp - w_poss)%w_max_w;
						if (buf_rev[buf_temp2].x % 100 < 80){ // sample randstrobes
							for (j = w_max_w+buf_temp2-w_min; j > w_max_w+buf_temp2-w_max; j--){
								uint64_t res = buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x;
								if (res < min_val){
									min_val = res;
									strobe_pos_next = j;
								}
							}
							info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
							info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max-w_poss)<<1 | 1;
							kv_push(mm128_t, km, *p, info_rev);
						} else { // sample k-mers
							// fprintf(stderr, "Rev: %d, %d\n", buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x, buf_rev[buf_temp2].x);
							info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x)) << 8 | strobe_len;
							info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max-w_poss)<<1 | 1;
							kv_push(mm128_t, km, *p, info_rev);
						}
					}
				}
				++w_buf_pos;
			}
		}
		if (++buf_pos == w_max_w) buf_pos = 0;
	}
    // compute rev minimizers on last w_max pos:
	mm128_t info_rev = { UINT64_MAX, UINT64_MAX };
	counter = 0;
	for (buf_temp = (w_max_w+buf_pos-w_max); buf_temp < (w_max_w+buf_pos); ++buf_temp){
		++counter;
		buf_temp2 = buf_temp%w_max_w;
		if (buf_rev[buf_temp2].x <= min.x) {
			min_val = UINT64_MAX;
			if (buf_rev[buf_temp2].x % 100 < 80){ // sample randstrobes
				for (j = w_max_w+buf_temp2-w_min; j > w_max_w+buf_temp2-w_max; j--){
					uint64_t res = buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x;
					if (res < min_val){
						min_val = res;
						strobe_pos_next = j;
					}
				}
				info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
				info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max+counter-1)<<1 | 1;
				kv_push(mm128_t, km, *p, info_rev);
			} else { // sample k-mers
				// fprintf(stderr, "Rev: %d, %d\n", buf_rev[(w_max_w+buf_pos-k)%w_max_w].x, buf_rev[buf_temp].x);
				info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x)) << 8 | strobe_len;
				info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max+counter-1)<<1 | 1;
				kv_push(mm128_t, km, *p, info_rev);
			}
			w_buf_pos = 0;
			min = buf_rev[buf_temp2];
		}
		if (++w_buf_pos == w){  // old min has moved outside the window
			// fprintf(stderr, "Moved outside window\n");
			min.x = UINT64_MAX;
			for (w_poss = (w-2); w_poss >= 0; --w_poss){
				// fprintf(stderr, "\tw_poss: %d, Fwd: %d, Rev: %d\n", w_poss, buf_fw[(buf_temp-w_poss)%w_max_w].x, buf_rev[(w_max_w+buf_pos-w_poss)%w_max_w].x);
				if (buf_rev[(buf_temp-w_poss)%w_max_w].x < min.x){
					w_buf_pos = w_poss;
					min = buf_rev[(buf_temp-w_poss)%w_max_w];
				}
			}
			for (w_poss = (w-2); w_poss >= 0; --w_poss){
				if (buf_rev[(buf_temp-w_poss)%w_max_w].x == min.x) {
					min_val = UINT64_MAX;
					buf_temp2 = (buf_temp - w_poss)%w_max_w;
					if (buf_rev[buf_temp2].x % 100 < 80){ // sample randstrobes
						for (j = w_max_w+buf_temp2-w_min; j > w_max_w+buf_temp2-w_max; j--){
							uint64_t res = buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x;
							if (res < min_val){
								min_val = res;
								strobe_pos_next = j;
							}
						}
						info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
						info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max-w_poss+counter-1)<<1 | 1;
						kv_push(mm128_t, km, *p, info_rev);
					} else { // sample k-mers
						// fprintf(stderr, "Rev: %d, %d\n", buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x, buf_rev[buf_temp2].x);
						info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[(w_max_w+buf_temp2-k)%w_max_w].x)) << 8 | strobe_len;
						info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max-w_poss+counter-1)<<1 | 1;
						kv_push(mm128_t, km, *p, info_rev);
					}
				}
			}
			++w_buf_pos;
		}
	}

}



/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param w_min  minimum strobe offset
 * @param w_max  maximum strobe offset
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch_randstrobes(void *km, const char *str, int len, int w, int k, int k_min, int w_min, int w_max, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, buf_temp, buf_temp2, min_pos, kmer_span, strobe_pos_next, w_buf_pos, w_poss, counter = 0;
	mm128_t buf_fw[256], min = { UINT64_MAX, UINT64_MAX };
	mm128_t buf_rev[256];
	tiny_queue_t tq;
	uint64_t min_val = UINT64_MAX;
	int w_max_w = 2*w_max + w;
	int strobe_len = w_max+k-1;
	int strobe_span = 0;
    int min_start = k+w_max+w;

	assert(len > 0 && (w_max_w > 0 && w_max_w < 256) && (k > 0 && k <= 14)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf_fw, 0xff, w_max_w * 16);
	memset(buf_rev, 0xff, w_max_w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info_fw = { UINT64_MAX, UINT64_MAX };
		mm128_t info_rev = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			// if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			// z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info_fw.x = hash64(kmer[0], mask); // << 8 | kmer_span;
				// info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
				info_rev.x = hash64(kmer[1], mask); // << 8 | kmer_span;
				// info_rev.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf_fw[buf_pos] = info_fw; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		buf_rev[buf_pos] = info_rev; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below

		if (l >= min_start){
			// fprintf(stderr, "Buf_pos: %d, w_buf %d, Min: %d, Fwd: %d, Rev: %d\n", buf_pos, w_buf_pos, min.x, buf_fw[(buf_pos+1+w)%w_max_w].x, info_rev.x);
			buf_temp = (w_max_w+buf_pos-w_max)%w_max_w;
			if ((buf_fw[buf_temp].x <= min.x) && (buf_fw[buf_temp].x <= buf_rev[buf_temp].x)){
				min_val = UINT64_MAX;
				info_fw.x = buf_fw[buf_temp].x;
				for (j = w_min+buf_temp; j < w_max+buf_temp; j++){
					uint64_t res = info_fw.x ^ buf_fw[j%w_max_w].x;
					if (res < min_val){
						min_val = res;
						strobe_pos_next = j;
					}
				}
				info_fw.x = ((info_fw.x << (2*k)) ^ (buf_fw[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
				info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
				kv_push(mm128_t, km, *p, info_fw);
				w_buf_pos = 0;
				min = buf_fw[buf_temp];
			}
			else if ((buf_rev[buf_temp].x <= min.x) && (l >= min_start+w_max)) {
				min_val = UINT64_MAX;
				for (j = w_max_w+buf_temp-w_min; j > w_max_w+buf_temp-w_max; j--){
					uint64_t res = buf_rev[buf_temp].x ^ buf_rev[j%w_max_w].x;
					if (res < min_val){
						min_val = res;
						strobe_pos_next = j;
					}
				}
				info_rev.x = ((buf_rev[buf_temp].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
				info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max)<<1 | 1;
				kv_push(mm128_t, km, *p, info_rev);
				w_buf_pos = 0;
				min = buf_rev[buf_temp];
			}

			if (++w_buf_pos == w){  // old min has moved outside the window
				// fprintf(stderr, "Moved outside window\n");
				min.x = UINT64_MAX;
                buf_temp = w_max_w+buf_pos-w_max;
				for (w_poss = (w-2); w_poss >= 0; --w_poss){
					// fprintf(stderr, "\tw_poss: %d, Fwd: %d, Rev: %d\n", w_poss, buf_fw[(buf_temp-w_poss)%w_max_w].x, buf_rev[(w_max_w+buf_pos-w_poss)%w_max_w].x);
					if (buf_fw[(buf_temp-w_poss)%w_max_w].x < min.x){
						w_buf_pos = w_poss;
						min = buf_fw[(buf_temp-w_poss)%w_max_w];
					}
					if ((buf_rev[(buf_temp-w_poss)%w_max_w].x < min.x) && (l >= min_start+w_max)){
						w_buf_pos = w_poss;
						min = buf_rev[(buf_temp-w_poss)%w_max_w];
					}
				}
				for (w_poss = (w-2); w_poss >= 0; --w_poss){
					if (buf_fw[(buf_temp-w_poss)%w_max_w].x == min.x){
						buf_temp2 = (buf_temp - w_poss)%w_max_w;
						min_val = UINT64_MAX;
						info_fw.x = buf_fw[buf_temp2].x;
						for (j = w_min+buf_temp2; j < w_max+buf_temp2; j++){
							uint64_t res = info_fw.x ^ buf_fw[j%w_max_w].x;
							if (res < min_val){
								min_val = res;
								strobe_pos_next = j;
							}
						}
						info_fw.x = ((info_fw.x << (2*k)) ^ (buf_fw[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
						info_fw.y = (uint64_t)rid<<32 | (uint32_t)(i-w_poss)<<1 | 0;
						kv_push(mm128_t, km, *p, info_fw);
					}
					if ((buf_rev[(buf_temp-w_poss)%w_max_w].x == min.x) && (l >= min_start+w_max)) {
						min_val = UINT64_MAX;
						buf_temp2 = (buf_temp - w_poss)%w_max_w;
						for (j = w_max_w+buf_temp2-w_min; j > w_max_w+buf_temp2-w_max; j--){
							uint64_t res = buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x;
							if (res < min_val){
								min_val = res;
								strobe_pos_next = j;
							}
						}
						info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
						info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max-w_poss)<<1 | 1;
						kv_push(mm128_t, km, *p, info_rev);
					}
				}
				++w_buf_pos;
			}
		}
		if (++buf_pos == w_max_w) buf_pos = 0;
	}
    // compute rev minimizers on last w_max pos:
	mm128_t info_rev = { UINT64_MAX, UINT64_MAX };
	counter = 0;
	for (buf_temp = (w_max_w+buf_pos-w_max); buf_temp < (w_max_w+buf_pos); ++buf_temp){
		++counter;
		buf_temp2 = buf_temp%w_max_w;
		if (buf_rev[buf_temp2].x <= min.x) {
			min_val = UINT64_MAX;
			for (j = w_max_w+buf_temp2-w_min; j > w_max_w+buf_temp2-w_max; j--){
				uint64_t res = buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x;
				if (res < min_val){
					min_val = res;
					strobe_pos_next = j;
				}
			}
			info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
			info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max+counter-1)<<1 | 1;
			kv_push(mm128_t, km, *p, info_rev);
			w_buf_pos = 0;
			min = buf_rev[buf_temp2];
		}
		if (++w_buf_pos == w){  // old min has moved outside the window
			// fprintf(stderr, "Moved outside window\n");
			min.x = UINT64_MAX;
			for (w_poss = (w-2); w_poss >= 0; --w_poss){
				// fprintf(stderr, "\tw_poss: %d, Fwd: %d, Rev: %d\n", w_poss, buf_fw[(buf_temp-w_poss)%w_max_w].x, buf_rev[(w_max_w+buf_pos-w_poss)%w_max_w].x);
				if (buf_rev[(buf_temp-w_poss)%w_max_w].x < min.x){
					w_buf_pos = w_poss;
					min = buf_rev[(buf_temp-w_poss)%w_max_w];
				}
			}
			for (w_poss = (w-2); w_poss >= 0; --w_poss){
				if (buf_rev[(buf_temp-w_poss)%w_max_w].x == min.x) {
					min_val = UINT64_MAX;
					buf_temp2 = (buf_temp - w_poss)%w_max_w;
					for (j = w_max_w+buf_temp2-w_min; j > w_max_w+buf_temp2-w_max; j--){
						uint64_t res = buf_rev[buf_temp2].x ^ buf_rev[j%w_max_w].x;
						if (res < min_val){
							min_val = res;
							strobe_pos_next = j;
						}
					}
					info_rev.x = ((buf_rev[buf_temp2].x << (2*k)) ^ (buf_rev[strobe_pos_next%w_max_w].x)) << 8 | strobe_len;
					info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-w_max-w_poss+counter-1)<<1 | 1;
					kv_push(mm128_t, km, *p, info_rev);
				}
			}
			++w_buf_pos;
		}
	}
}


u_int64_t seq_to_hashvalue_fw(const char *str, int start, int length)
{
	uint64_t kmask = (1ULL<<2*length) - 1;
	uint64_t x = 0;
    for (int i = start; i < start+length; i++) {
        int c = seq_nt4_table[(uint8_t) str[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;
            }
        else {
            return 0;
        }
    }
    uint64_t hash_k = hash64(x, kmask);
    return hash_k;
}


u_int64_t seq_to_hashvalue_rev(const char *str, int start, int length)
{
	uint64_t kmask = (1ULL<<2*length) - 1;
	uint64_t kshift = 2 * (length - 1);
	uint64_t x = 0;
    for (int i = start-length; i < start; i++) {
        int c = seq_nt4_table[(uint8_t) str[i]];
        if (c < 4) { // not an "N" base
			x = (x >> 2) | (3ULL^c) << kshift;
            }
        else {
            return 0;
        }
    }
    uint64_t hash_k = hash64(x, kmask);
    return hash_k;
}


/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param k_min  minimum k-mer size
 * @param w_min  minimum strobe offset
 * @param w_max  maximum strobe offset
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch_multistrobes(void *km, const char *str, int len, int w, int k, int k_min, int w_min, int w_max, uint32_t rid, int is_hpc, mm128_v *p)
{
	// fprintf(stderr, "k: %d, k_min: %d, w_min: %d, w_max: %d\n\n", 2*k, k_min, w_min, w_max);
	uint64_t shift1 = 2 * (k_min - 1), mask = (1ULL<<2*k_min) - 1, kmer[2] = {0,0};
	int i, j, l, k1, offset, buf_pos, buf_temp, buf_temp2, min_pos, kmer_span, strobe_pos_next, w_buf_pos, w_poss, counter = 0;
	uint64_t hash1 = 0, hash2 = 0;
	mm128_t buf_fw[256], min = { UINT64_MAX, UINT64_MAX };
	mm128_t buf_rev[256];
	tiny_queue_t tq;
	uint64_t min_val = UINT64_MAX;
	int w_max_w = 2*w_max + 2*(2*k-k_min)+w;
	int strobe_len = w_max + k - k_min-1;
	int strobe_span = 0;
    int min_start = k+w_max+w;
	int shift = w_max + k - k_min;

	int w_min_off = w_min - k; // + k1
	int w_max_off = w_max - k; // + k1

	assert(len > 0 && (w_max_w > 0 && w_max_w < 256) && (k_min > 0) && (k > 0 && k <= 14)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf_fw, 0xff, w_max_w * 16);
	memset(buf_rev, 0xff, w_max_w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	k = 2*k;
	int nr_options = k-2*k_min+1;

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info_fw = { UINT64_MAX, UINT64_MAX };
		mm128_t info_rev = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k_min) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k_min? l + 1 : k_min;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			// if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			// z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k_min && kmer_span < 256) {
				info_fw.x = hash64(kmer[0], mask); // << 8 | kmer_span;
				// info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
				info_rev.x = hash64(kmer[1], mask); // << 8 | kmer_span;
				// info_rev.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf_fw[buf_pos] = info_fw; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		buf_rev[buf_pos] = info_rev; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below

		if (l >= min_start){
			// fprintf(stderr, "Buf_pos: %d, w_buf %d, Min: %d, Fwd: %d, Rev: %d\n", buf_pos, w_buf_pos, min.x, buf_fw[(buf_pos+1+w)%w_max_w].x, info_rev.x);
			buf_temp = (w_max_w+buf_pos-shift)%w_max_w;
			if ((buf_fw[buf_temp].x <= min.x) && (buf_fw[buf_temp].x <= buf_rev[buf_temp].x)){
				k1 = k_min + buf_fw[buf_temp].x % nr_options;
				min_val = UINT64_MAX;
				hash1 = seq_to_hashvalue_fw(str, i-shift-k_min+1, k1);
				int offset = w_min_off+k1;
				for (j = w_min_off+buf_temp+k1; j < w_max_off+buf_temp+k1; j++){
					uint64_t res = hash1 ^ buf_fw[j%w_max_w].x;
					if (res < min_val){
						min_val = res;
						// strobe_pos_next = j;
						hash2 = offset;
					}
					++offset;
				}
				hash2 = seq_to_hashvalue_fw(str, hash2+i-shift-k_min+1, k-k1);
				info_fw.x = ((hash1 << (2*k1)) ^ hash2) << 8 | strobe_len;
				info_fw.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
				kv_push(mm128_t, km, *p, info_fw);
				w_buf_pos = 0;
				min = buf_fw[buf_temp];
			}
			else if ((buf_rev[buf_temp].x <= min.x) && (l >= min_start+w_max)) {
				k1 = k_min + buf_rev[buf_temp].x % nr_options;
				hash1 = seq_to_hashvalue_rev(str, i-shift+1, k1);
				min_val = UINT64_MAX;
				int offset = w_max_w-w_min_off-k1;
				for (j = w_max_w+buf_temp-w_min_off-k1; j > w_max_w+buf_temp-w_max_off-k1; j--){
					uint64_t res = hash1 ^ buf_rev[j%w_max_w].x;
					if (res < min_val){
						min_val = res;
						// strobe_pos_next = j;
						hash2 = offset;
					}
					--offset;
				}
				hash2 = seq_to_hashvalue_rev(str, hash2-w_max_w+i-shift+1, k-k1);
				info_rev.x = ((hash1 << (2*k1)) ^ hash2) << 8 | strobe_len;
				info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-shift)<<1 | 1;
				kv_push(mm128_t, km, *p, info_rev);
				w_buf_pos = 0;
				min = buf_rev[buf_temp];
			}

			if (++w_buf_pos == w){  // old min has moved outside the window
				// fprintf(stderr, "Moved outside window\n");
				min.x = UINT64_MAX;
                buf_temp = w_max_w+buf_pos-shift;
				for (w_poss = (w-2); w_poss >= 0; --w_poss){
					// fprintf(stderr, "\tw_poss: %d, Fwd: %d, Rev: %d\n", w_poss, buf_fw[(buf_temp-w_poss)%w_max_w].x, buf_rev[(w_max_w+buf_pos-w_poss)%w_max_w].x);
					if (buf_fw[(buf_temp-w_poss)%w_max_w].x < min.x){
						w_buf_pos = w_poss;
						min = buf_fw[(buf_temp-w_poss)%w_max_w];
					}
					if ((buf_rev[(buf_temp-w_poss)%w_max_w].x < min.x) && (l >= min_start+w_max+k_min)){
						w_buf_pos = w_poss;
						min = buf_rev[(buf_temp-w_poss)%w_max_w];
					}
				}
				for (w_poss = (w-2); w_poss >= 0; --w_poss){
					if (buf_fw[(buf_temp-w_poss)%w_max_w].x == min.x){
						buf_temp2 = (buf_temp - w_poss)%w_max_w;
						k1 = k_min + buf_fw[buf_temp2].x % nr_options;
						hash1 = seq_to_hashvalue_fw(str, i-w_poss-shift-k_min+1, k1);
						min_val = UINT64_MAX;
						int offset = w_min_off+k1;
						for (j = w_min_off+buf_temp2+k1; j < w_max_off+buf_temp2+k1; j++){
							uint64_t res = hash1 ^ buf_fw[j%w_max_w].x;
							if (res < min_val){
								min_val = res;
								// strobe_pos_next = j;
								hash2 = offset;
							}
							++offset;
						}
						hash2 = seq_to_hashvalue_fw(str, hash2+i-w_poss-shift-k_min+1, k-k1);
						info_fw.x = ((hash1 << (2*k1)) ^ hash2) << 8 | strobe_len;
						info_fw.y = (uint64_t)rid<<32 | (uint32_t)(i-w_poss)<<1 | 0;
						kv_push(mm128_t, km, *p, info_fw);
					}
					if ((buf_rev[(buf_temp-w_poss)%w_max_w].x == min.x) && (l >= min_start+w_max+k_min)) {
						buf_temp2 = (buf_temp - w_poss)%w_max_w;
						k1 = k_min + buf_rev[buf_temp2].x % nr_options;
						hash1 = seq_to_hashvalue_rev(str, i-w_poss-shift+1, k1);
						min_val = UINT64_MAX;
						int offset = w_max_w-w_min_off-k1;
						for (j = w_max_w+buf_temp2-w_min_off-k1; j > w_max_w+buf_temp2-w_max_off-k1; j--){
							uint64_t res = hash1 ^ buf_rev[j%w_max_w].x;
							if (res < min_val){
								min_val = res;
								// strobe_pos_next = j;
								hash2 = offset;
							}
							--offset;
						}
						hash2 = seq_to_hashvalue_rev(str, hash2-w_max_w+i-w_poss-shift+1, k-k1);
						info_rev.x = ((hash1 << (2*k1)) ^ hash2) << 8 | strobe_len;
						info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-shift-w_poss)<<1 | 1;
						kv_push(mm128_t, km, *p, info_rev);
					}
				}
				++w_buf_pos;
			}
		}
		if (++buf_pos == w_max_w) buf_pos = 0;
	}
    // compute rev minimizers on last w_max pos:
	mm128_t info_rev = { UINT64_MAX, UINT64_MAX };
	counter = 0;
	for (buf_temp = (w_max_w+buf_pos-shift); buf_temp < (w_max_w+buf_pos); ++buf_temp){
		++counter;
		buf_temp2 = buf_temp%w_max_w;
		if (buf_rev[buf_temp2].x <= min.x) {
			k1 = k_min + buf_rev[buf_temp2].x % nr_options;
			hash1 = seq_to_hashvalue_rev(str, i+counter-shift, k1);
			min_val = UINT64_MAX;
			int offset = w_max_w-w_min_off-k1;
			for (j = w_max_w+buf_temp2-w_min_off-k1; j > w_max_w+buf_temp2-w_max_off-k1; j--){
				uint64_t res = hash1 ^ buf_rev[j%w_max_w].x;
				if (res < min_val){
					min_val = res;
					// strobe_pos_next = j;
					hash2 = offset;
				}
				--offset;
			}
			hash2 = seq_to_hashvalue_rev(str, hash2-w_max_w+i+counter-shift, k-k1);
			info_rev.x = ((hash1 << (2*k1)) ^ hash2) << 8 | strobe_len;
			info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-shift+counter-1)<<1 | 1;
			kv_push(mm128_t, km, *p, info_rev);
			w_buf_pos = 0;
			min = buf_rev[buf_temp2];
		}
		if (++w_buf_pos == w){  // old min has moved outside the window
			// fprintf(stderr, "Moved outside window\n");
			min.x = UINT64_MAX;
			for (w_poss = (w-2); w_poss >= 0; --w_poss){
				// fprintf(stderr, "\tw_poss: %d, Fwd: %d, Rev: %d\n", w_poss, buf_fw[(buf_temp-w_poss)%w_max_w].x, buf_rev[(w_max_w+buf_pos-w_poss)%w_max_w].x);
				if (buf_rev[(buf_temp-w_poss)%w_max_w].x < min.x){
					w_buf_pos = w_poss;
					min = buf_rev[(buf_temp-w_poss)%w_max_w];
				}
			}
			for (w_poss = (w-2); w_poss >= 0; --w_poss){
				if (buf_rev[(buf_temp-w_poss)%w_max_w].x == min.x) {
					buf_temp2 = (buf_temp - w_poss)%w_max_w;
					k1 = k_min + buf_rev[buf_temp2].x % nr_options;
					hash1 = seq_to_hashvalue_rev(str, i+counter-w_poss-shift, k1);
					min_val = UINT64_MAX;
					int offset = w_max_w-w_min_off-k1;
					for (j = w_max_w+buf_temp2-w_min_off-k1; j > w_max_w+buf_temp2-w_max_off-k1; j--){
						uint64_t res = hash1 ^ buf_rev[j%w_max_w].x;
						if (res < min_val){
							min_val = res;
							// strobe_pos_next = j;
							hash2 = offset;
						}
						--offset;
					}
					hash2 = seq_to_hashvalue_rev(str, hash2-w_max_w+i+counter-w_poss-shift, k-k1);
					info_rev.x = ((hash1 << (2*k1)) ^ hash2) << 8 | strobe_len;
					info_rev.y = (uint64_t)rid<<32 | (uint32_t)(i-shift-w_poss+counter-1)<<1 | 1;
					kv_push(mm128_t, km, *p, info_rev);
				}
			}
			++w_buf_pos;
		}
	}
}