/*
 * 
 * Copyright (c) 2013, Zechen Chong <chongzechen@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FILTER_KMER_H
#define FILTER_KMER_H

#include <stdint.h>
#include <unistd.h>
#include "timer.h"
#include "file_reader.h"
#include "hashset.h"
#include "string.h"
#include "list.h"
#include "dna.h"
#include "sort.h"
#include "bitvec.h"
#include "counting_bloom_filter.h"

#define UNIQ_KMER_MAX_KSIZE 31
#define UNIQ_KMER_MAX_CNT 999

typedef struct {
	uint64_t kmer, cnt:16, cnt2:16;
} kmer_t;

typedef struct {
	char *name;
	char *header;
	char *seq;
	char *qual;
} pair;

typedef enum {SOMATIC, GERMLINE, LOH} mut_type;

typedef struct {
	uint32_t pid;
	mut_type mt;

} pair_t;

define_list(pairv, pair_t);
define_list(flist, char*);


#define kmer_hashcode(k) u64hashcode((k).kmer)
#define kmer_equals(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(kmerhash, kmer_t, kmer_hashcode, kmer_equals);

static inline int trim_lowq(char *qual, int len, char min) {
	int i;
	for (i = 0; i < len; i++) {
		if (qual[i] <= min) {
			if (i+1 < len && qual[i+1] <= min) {
				break;
			}
		}
	}
	return i;
}

#ifdef __CPLUSPLUS
extern "C" {
#endif

BitVec* build_kmerhash(FileReader *fr, uint32_t ksize, int is_fq, BitVec *bt, u64hash *refhash, uint64_t *idx);
kmerhash* build_readshash(FileReader *readfr, uint32_t ksize, int is_fq, kmerhash *hash, BitVec *bt, uint64_t *idx, CBF *occ_table, uint32_t mincnt);
kmerhash* build_refkmerhash(FileReader *fr, FileReader *read1fr, FileReader *read2fr, int is_fq, uint32_t ksize, kmerhash* hash, uint32_t mincnt);
void cal_ctrl_kmers(kmerhash *hash, FileReader *fr, uint32_t ksize, int is_fq);
uint32_t count_readnum(FileReader *readfr, int is_fq);
uint64_t filter_ref_kmers(kmerhash *hash, FileReader *fr, uint32_t ksize);
pairv* loadkmerseq(kmerhash *hash, uint32_t ksize, uint32_t mincnt, uint32_t maxcnt2, FileReader *f1, FileReader *f2, int is_somatic);
void dedup_pairs(pairv *pairs, FILE *out1, FILE *out2, FILE *out3, FILE *out4, FileReader *f1, FileReader *f2);
void destroy_pairv(pairv *pairs);

#ifdef __CPLUSPLUS
}
#endif

#endif
