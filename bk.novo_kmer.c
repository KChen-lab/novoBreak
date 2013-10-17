/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
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

#include "string.h"
#include "file_reader.h"
#include "list.h"
#include "dna.h"
#include "sort.h"
#include "bitvec.h"
#include "hashset.h"
#include <stdint.h>
#include <unistd.h>

#define UNIQ_KMER_MAX_KSIZE 27
#define UNIQ_KMER_MAX_CNT 1023

typedef struct {
	uint64_t kmer:54, cnt:10;
} Kmer;

#define kmer_hashcode(k) u64hashcode((k).kmer)
#define kmer_equals(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(kmerhash, Kmer, kmer_hashcode, kmer_equals);

define_list(strv, String*);
define_list(u32vv, u32list*);
define_list(u64vv, u64list*);

typedef struct {
	strv    *names;
	u64vv   *seqs;
	u32list *seqlens;
	u32vv   *muts;
	cuhash  *index;
} RefSeq;

typedef struct {
	uint64_t seqoff:44;
	uint64_t seqlen1:8, seqlen2:8;
	uint64_t chr:16, dir:1, dir1:1, dir2:1;
	uint64_t novo:1, lr:1, pos:32, off:15;
} PairHit;

define_list_core(pairv, PairHit, uint32_t, 0xFU);

typedef struct {
	uint32_t chr1:16, chr2:16;
	uint32_t pos1:31, dir1:1, pos2:31, dir2:1;
} MateHit;

define_list_core(matev, MateHit, uint32_t, 0xFU);

typedef struct {
	uint64_t kmer:54, cnt:10;
	pairv *pairs;
	matev *mates;
} NovoKmer;

typedef struct {
	uint32_t id, chr:16, off, len:15, dir:1;
} Range;

define_list_core(rangev, Range, uint32_t, 0xFU);

typedef struct {
	uint32_t bin_size;
	BitVec **bins;
	rangev ***ranges;
} MateBin;

define_list(nkv, NovoKmer);

#define novo_hashcode(k) u64hashcode((k).kmer)
#define novo_equals(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(novohash, NovoKmer, novo_hashcode, novo_equals);

int check_sam_cigar(char *cigar, int len){
	int i;
	for(i=0;i<len;i++){
		if(!((cigar[i] >= '0' && cigar[i] <= '9') || cigar[i] == 'M')) return 0;
	}
	return 1;
}

kmerhash* build_kmerhash(FileReader *fr, uint32_t ksize){
	kmerhash *hash;
	Kmer KMER, *kmer;
	char *seq, *str;
	uint64_t k, r, kmask;
	uint32_t i, len, f, flag, rid;
	int exists, m, n, nm, skip;
	hash = init_kmerhash(1023);
	KMER.cnt = 0;
	KMER.kmer = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	flag = 0;
	rid = 0;
	while((n = fread_table(fr)) != -1){
		if(n < 11) continue;
		if(fr->line->string[0] == '#') continue;
		// checking SAM file
		f = atoi(get_col_str(fr, 1));
		if(flag == 0){
			if((f & 0x40) == 0) flag = 3;
			else flag = 1;
		} else if(flag == 1){
			if((f & 0x80) == 0) flag = 3;
			else flag = 0;
		}
		if(flag == 3){
			fprintf(stderr, " Two lines in SAM file are not properly paired, please check your SAM file -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
			fprintf(stderr, "%s\n", fr->line->string);
			fflush(stderr);
			abort();
		}
		if((rid & 0xFFFFU) == 0){
			fprintf(stdout, "[%s] parsed %10u line\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		rid ++;
		skip = 0;
		if(check_sam_cigar(get_col_str(fr, 5), get_col_len(fr, 5))){
			for(m=11;m<n;m++){
				if(get_col_len(fr, m) < 6) continue;
				str = get_col_str(fr, m);
				if(memcmp(str, "NM:i:", 6) != 0) continue;
				nm = atoi(str + 6);
				if(nm <= 1) skip = 1;
				break;
			}
		}
		if(skip) continue;
		seq = get_col_str(fr, 9);
		len = get_col_len(fr, 9);
		k = 0;
		for(i=0;i<len;i++){
			k = ((k << 2) | base_bit_table[(int)seq[i]]) & kmask;
			if(i + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if(r < k){
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			kmer = prepare_kmerhash(hash, KMER, &exists);
			if(exists){
				if(kmer->cnt < UNIQ_KMER_MAX_CNT) kmer->cnt ++;
			} else {
				kmer->kmer = KMER.kmer;
				kmer->cnt  = 1;
			}
		}
	}
	fprintf(stdout, "[%s] parsed %10u line\n", __FUNCTION__, rid);
	fflush(stdout);
	return hash;
}

RefSeq* build_refseq(FileReader *ref, FileReader *mut, uint32_t ksize){
	RefSeq *rs;
	String *str;
	Sequence *seq;
	u32list *muts;
	u64list *seqs;
	cuhash_t CU, *cu;
	char chs[256];
	uint32_t v, vv, len;
	int i, n;
	rs = malloc(sizeof(RefSeq));
	rs->names   = init_strv(1024);
	rs->seqs    = init_u64vv(1024);
	rs->seqlens = init_u32list(1024);
	rs->muts    = init_u32vv(1024);
	rs->index   = init_cuhash(1023);
	seq = NULL;
	while(fread_fasta(&seq, ref)){
		str = clone_string(&(seq->name));
		CU.key = str->string;
		CU.val = count_strv(rs->names);
		put_cuhash(rs->index, CU);
		push_strv(rs->names, str);
		seqs = init_u64list(seq->seq.size / 32 + 2);
		seq2bits(ref_u64list(seqs, 0), 0, seq->seq.string, seq->seq.size);
		push_u64vv(rs->seqs, seqs);
		push_u32list(rs->seqlens, seq->seq.size);
		muts = init_u32list(1024);
		push_u32vv(rs->muts, muts);
	}
	if(mut == NULL) return rs;
	//chr1\t10222\tC\tG\tA
	muts = NULL;
	len = 0;
	chs[0] = 0;
	while((n = fread_table(mut)) != -1){
		if(n < 3) continue;
		if(strcmp(chs, get_col_str(mut, 0)) != 0){
			strcpy(chs, get_col_str(mut, 0));
			CU.key = chs;
			CU.val = 0;
			cu = get_cuhash(rs->index, CU);
			if(cu == NULL){
				muts = NULL;
				fprintf(stderr, " -- Ignore chromosome[%s] in %s -- %s:%d --\n", chs, __FUNCTION__, __FILE__, __LINE__);
				fflush(stderr);
				len = 0;
			} else {
				muts = get_u32vv(rs->muts, cu->val);
				len = get_u32list(rs->seqlens, cu->val);
			}
		}
		if(muts == NULL) continue;
		v = atoi(get_col_str(mut, 1));
		if(v < ksize || v + ksize > len) continue;
		v <<= 4;
		for(i=2;i<n;i++){
			vv = base_bit_table[(int)get_col_str(mut, i)[0]];
			if(vv == 4) continue;
			v |= (1U << vv);
		}
		push_u32list(muts, v);
	}
	return rs;
}

void free_refseq(RefSeq *rs){
	uint32_t i;
	for(i=0;i<count_strv(rs->names);i++) free_string(get_strv(rs->names, i));
	for(i=0;i<count_u64vv(rs->seqs);i++) free_u64list(get_u64vv(rs->seqs, i));
	free_u32list(rs->seqlens);
	for(i=0;i<count_u32vv(rs->muts);i++) free_u32list(get_u32vv(rs->muts, i));
	free_cuhash(rs->index);
	free(rs);
}

uint64_t filter_ref_genome(kmerhash *hash, RefSeq *rs, uint32_t ksize){
	Kmer KMER;
	uint64_t *seqs;
	uint64_t k, r, kmask, ret;
	uint32_t i, j, len;
	KMER.cnt = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	ret = 0;
	for(j=0;j<count_u64vv(rs->seqs);j++){
		seqs = ref_u64list(get_u64vv(rs->seqs, j), 0);
		len  = get_u32list(rs->seqlens, j);
		k = 0;
		fprintf(stdout, "[%s] %s\n", __FUNCTION__, get_strv(rs->names, j)->string);
		fflush(stdout);
		for(i=0;i<len;i++){
			k = ((k << 2) | bits2bit(seqs, i)) & kmask;
			if(i + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if(r < k){
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			ret += remove_kmerhash(hash, KMER);
		}
	}
	return ret;
}

uint64_t filter_separate_dbSNP(kmerhash *hash, RefSeq *rs, uint32_t ksize){
	Kmer KMER;
	uint64_t *seqs;
	u32list *muts;
	uint64_t k, r, kmask, ret;
	uint32_t i, j, s, e, len, v, f;
	KMER.cnt = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	ret = 0;
	for(i=0;i<count_u32vv(rs->muts);i++){
		muts = get_u32vv(rs->muts, i);
		seqs = ref_u64list(get_u64vv(rs->seqs, i), 0);
		len  = get_u32list(rs->seqlens, i);
		fprintf(stdout, "[%s] %s\n", __FUNCTION__, get_strv(rs->names, i)->string);
		fflush(stdout);
		for(j=0;j<count_u32list(muts);j++){
			v = get_u32list(muts, j);
			f = v & 0x0FU;
			v = v >> 4;
			//if(v < ksize && v + ksize > len) continue;
			v = v - ksize;
			for(e=0;e<4;e++){
				if(((f >> e) & 0x01) == 0) continue;
				k = 0;
				for(s=0;s+1<ksize;s++){ k = (k << 2) | bits2bit(seqs, v - ksize + s); }
				k = (k << 2) | e;
				for(s=1;s<=ksize;s++){
					r = dna_rev_seq(k, ksize);
					KMER.kmer = (r < k)? r : k;
					ret += remove_kmerhash(hash, KMER);
					k = ((k << 2) | bits2bit(seqs, v + s)) & kmask;
				}
			}
		}
	}
	return ret;
}

uint64_t filter_combine_dbSNP(kmerhash *hash, RefSeq *rs, uint32_t ksize){
	Kmer KMER;
	uint64_t *seqs;
	u32list *muts;
	uint64_t kk, k, r, kmask, ret;
	uint32_t i, j, s, t, len;
	uint32_t x, y, a, b, es[32], fs[32];
	KMER.cnt = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	ret = 0;
	for(i=0;i<count_u32vv(rs->muts);i++){
		muts = get_u32vv(rs->muts, i);
		seqs = ref_u64list(get_u64vv(rs->seqs, i), 0);
		len  = get_u32list(rs->seqlens, i);
		fprintf(stdout, "[%s] %s\n", __FUNCTION__, get_strv(rs->names, i)->string);
		fflush(stdout);
		for(x=0,y=0;x+1<count_u32list(muts);x++){
			if(y < x) y = x;
			while(y + 1 < count_u32list(muts) && (get_u32list(muts, y + 1) >> 4) < (get_u32list(muts, x) >> 4) + ksize) y ++;
			if(y == x) continue;
			for(j=x;j<=y;j++) es[j-x] = get_u32list(muts, j);
			for(j=2;j<=y-x+1;j++){
				a = es[j - 1] >> 4;
				b = (j == y - x + 1)? (es[0] >> 4) + ksize : (es[j] >> 4);
				kk = 0;
				for(s=0;s+1<ksize;s++) kk = (kk << 2) | bits2bit(seqs, a + s - ksize);
				for(s=a;s<b;s++){
					kk = ((kk << 2) | bits2bit(seqs, s)) & kmask;
					for(t=0;t<j;t++) fs[t] = 0;
					while(1){
						t = j;
						while(t){
							do { fs[t-1] ++; } while(fs[t-1] < 4 && ((es[t-1] >> fs[t-1]) & 0x01) == 0);
							if(t > 1 && fs[t-1] == 4){
								fs[t-1] = 0;
								while(fs[t-1] < 4 && ((es[t-1] >> fs[t-1]) & 0x01) == 0){ fs[t-1] ++; }
								t --;
							} else break;
						}
						//for(t=0;t<j;t++) fprintf(stdout, "[%u][%u] %u\n", s, t, fs[t]);
						if(fs[0] >= 4) break;
						k = kk;
						for(t=0;t<j;t++){
							k &= ~(0x03LLU << (((es[t] >> 4) - s) << 1));
							k |= fs[t] << (((es[t] >> 4) - s) << 1);
						}
						r = dna_rev_seq(k, ksize);
						KMER.kmer = (r < k)? r : k;
						ret += remove_kmerhash(hash, KMER);
					}
				}
			}
		}
	}
	return ret;
}

uint64_t filter_control_kmers(kmerhash *hash, FileReader *fr, uint32_t ksize){
	Kmer KMER;
	char *seq, *str;
	uint64_t k, r, kmask, ret;
	uint32_t i, len, rid;
	int n, m, nm, skip;
	KMER.cnt = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	ret = 0;
	rid = 0;
	while((n = fread_table(fr)) != -1){
		if(n < 11) continue;
		if(fr->line->string[0] == '#') continue;
		if((rid & 0xFFFFU) == 0){
			fprintf(stdout, "[%s] parsed %10u line\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		rid ++;
		skip = 0;
		if(check_sam_cigar(get_col_str(fr, 5), get_col_len(fr, 5))){
			for(m=11;m<n;m++){
				if(get_col_len(fr, m) < 6) continue;
				str = get_col_str(fr, m);
				if(memcmp(str, "NM:i:", 6) != 0) continue;
				nm = atoi(str + 6);
				if(nm <= 1)skip = 1;
				break;
			}
		}
		if(skip) continue;
		seq = get_col_str(fr, 9);
		len = get_col_len(fr, 9);
		k = 0;
		for(i=0;i<len;i++){
			k = ((k << 2) | base_bit_table[(int)seq[i]]) & kmask;
			if(i + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if(r < k){
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			ret += remove_kmerhash(hash, KMER);
		}
	}
	fprintf(stdout, "[%s] parsed %10u line\n", __FUNCTION__, rid);
	fflush(stdout);
	return ret;
}

novohash* kmerhash2novohash(kmerhash *khash, uint32_t min_cnt, uint32_t max_cnt){
	novohash *nhash;
	Kmer *kmer;
	NovoKmer NK;
	nhash = init_novohash(1023);
	while((kmer = ref_iter_kmerhash(khash))){
		if(kmer->cnt < min_cnt || kmer->cnt > max_cnt) continue;
		NK.kmer  = kmer->kmer;
		NK.cnt   = kmer->cnt;
		NK.pairs = init_pairv(4);
		NK.mates = init_matev(4);
		put_novohash(nhash, NK);
	}
	free_kmerhash(khash);
	return nhash;
}

int cmp_pairhit(PairHit h1, PairHit h2, void *obj){
	if(h1.lr == h2.lr){
		if(h1.chr == h2.chr){
			return ((int)h1.pos) - ((int)h2.pos);
		} else if(h1.chr < h2.chr) return -1;
		else return 1;
	} else if(h1.lr) return -1;
	else return 1;
	obj = obj;
}

define_quick_sort(sort_pairhits, PairHit, cmp_pairhit);

#define cmp_range_by_id(r1, r2, obj) (int)(((int64_t)(r1).id) - ((int64_t)(r2).id))

define_quick_sort(sort_range_by_id, Range, cmp_range_by_id);

#define get_tabs_str(tabs, idx) ((VirtualString*)get_vec_ref(tabs, idx))->string
#define get_tabs_len(tabs, idx) ((VirtualString*)get_vec_ref(tabs, idx))->size

inline uint32_t abs_sub(uint32_t a, uint32_t b){ return (a < b)? b - a : a - b; }

uint64_t load_novo_kmer_pairs(novohash *hash, RefSeq *rs, FileReader *fr, uint32_t ksize, u64list *bits, uint32_t *avg_ins){
	NovoKmer NK, *nk;
	PairHit *pair;
	Vector *tabs1, *tabs2;
	String *line1, *line2;
	VirtualString vl1, vl2;
	cuhash_t CU, *cu;
	uint64_t ret, offset, k, r, kmask, t_ins, n_ins;
	uint32_t i, len, rid, rmask, used, flag1, flag2, uniq1, uniq2, chr1, chr2;
	char *seq, *str;
	int m, n, nm, skip;
	ret = 0;
	NK.pairs = NULL;
	NK.mates = NULL;
	NK.cnt = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	offset = 0;
	rid = 0;
	line1 = init_string(1024);
	line2 = init_string(1024);
	tabs1 = init_vec(sizeof(VirtualString), 32);
	tabs2 = init_vec(sizeof(VirtualString), 32);
	t_ins = 0;
	n_ins = 0;
	while(1){
		chr1 = chr2 = 0xFFFFU;
		while((n = fread_line(line1, fr)) != -1){
			if(line1->string[0] == '#') continue;
			vl1.string = line1->string;
			vl1.size   = line1->size;
			clear_vec(tabs1);
			n = split_vstring(&vl1, '\t', tabs1, 1);
			if(n < 11) continue;
			flag1 = atoi(get_tabs_str(tabs1, 1));
			uniq1 = 1;
			if(flag1 & 0x04) uniq1 = 0;
			else {
				for(m=11;m<n;m++){
					if(get_tabs_len(tabs1, m) < 6) continue;
					if(memcmp(get_tabs_str(tabs1, m), "XT:A:R", 6) == 0){ uniq1 = 2; break; }
				}
			}
			if(uniq1 == 1){
				CU.key = get_tabs_str(tabs1, 2);
				CU.val = 0;
				cu = get_cuhash(rs->index, CU);
				if(cu == NULL) chr1 = 0xFFFFU;
				else chr1 = cu->val;
			} else chr1 = 0xFFFFU;
			rid ++;
			break;
		}
		if(n == -1) break;
		while((n = fread_line(line2, fr)) != -1){
			if(line2->string[0] == '#') continue;
			vl2.string = line2->string;
			vl2.size   = line2->size;
			clear_vec(tabs2);
			n = split_vstring(&vl2, '\t', tabs2, 1);
			if(n < 11) continue;
			flag2 = atoi(get_tabs_str(tabs2, 1));
			uniq2 = 1;
			if(flag2 & 0x04) uniq2 = 0;
			else {
				for(m=11;m<n;m++){
					if(get_tabs_len(tabs2, m) < 6) continue;
					if(memcmp(get_tabs_str(tabs2, m), "XT:A:R", 6) == 0){ uniq2 = 2; break; }
				}
			}
			if(uniq2 == 1){
				CU.key = get_tabs_str(tabs2, 2);
				CU.val = 0;
				cu = get_cuhash(rs->index, CU);
				if(cu == NULL) chr2 = 0xFFFFU;
				else chr2 = cu->val;
			} else chr2 = 0xFFFFU;
			rid ++;
			break;
		}
		if(n == -1) break;
		if((rid & 0xFFFFU) == 0){
			fprintf(stdout, "[%s] parsed %10u line\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		n = abs_sub(atol(get_tabs_str(tabs1, 3)), atol(get_tabs_str(tabs2, 3)));
		if(chr1 == chr2 && chr1 != 0xFFFFU && n < 1024 * 1024){
			t_ins += n;
			n_ins ++;
		}
		used = 0;
		for(n=0;n<1;n++){
			skip = 0;
			if(check_sam_cigar(get_tabs_str(tabs1, 5), get_tabs_len(tabs1, 5))){
				for(m=11;m<(int)vec_size(tabs1);m++){
					if(get_tabs_len(tabs1, m) < 6) continue;
					str = get_tabs_str(tabs1, m);
					if(memcmp(str, "NM:i:", 6) != 0) continue;
					nm = atoi(str + 6);
					if(nm <= 1) skip = 1;
					break;
				}
			}
			if(skip) continue;
			seq = get_tabs_str(tabs1, 9);
			len = get_tabs_len(tabs1, 9);
			k = 0;
			for(i=0;i<len;i++){
				k = ((k << 2) | base_bit_table[(int)seq[i]]) & kmask;
				if(i + 1 < ksize) continue;
				r = dna_rev_seq(k, ksize);
				if(r < k){
					NK.kmer = r;
					rmask = 1;
				} else {
					NK.kmer = k;
					rmask = 0;
				}
				nk = get_novohash(hash, NK);
				if(nk == NULL) continue;
				pair = next_ref_pairv(nk->pairs);
				pair->seqoff  = offset;
				pair->seqlen1 = get_tabs_len(tabs1, 9);
				pair->seqlen2 = get_tabs_len(tabs2, 9);
				pair->chr     = chr2;
				pair->lr      = rmask ^ ((flag1 & 0x10)? 1:0);
				pair->dir1    = rmask;
				pair->dir2    = 1 ^ pair->lr ^ ((flag2 & 0x10)? 1:0);
				pair->dir     = pair->dir2;
				pair->novo    = 0;
				pair->pos     = atol(get_tabs_str(tabs2, 3));
				pair->off     = rmask? len - i - 1 : i + 1 - ksize;
				used          = 1;
			}
		}
		for(n=0;n<1;n++){
			skip = 0;
			if(check_sam_cigar(get_tabs_str(tabs2, 5), get_tabs_len(tabs2, 5))){
				for(m=11;m<(int)vec_size(tabs2);m++){
					if(get_tabs_len(tabs2, m) < 6) continue;
					str = get_tabs_str(tabs2, m);
					if(memcmp(str, "NM:i:", 6) != 0) continue;
					nm = atoi(str + 6);
					if(nm <= 1) skip = 1;
					break;
				}
			}
			if(skip) continue;
			seq = get_tabs_str(tabs2, 9);
			len = get_tabs_len(tabs2, 9);
			k = 0;
			for(i=0;i<len;i++){
				k = ((k << 2) | base_bit_table[(int)seq[i]]) & kmask;
				if(i + 1 < ksize) continue;
				r = dna_rev_seq(k, ksize);
				if(r < k){
					NK.kmer = r;
					rmask = 1;
				} else {
					NK.kmer = k;
					rmask = 0;
				}
				nk = get_novohash(hash, NK);
				if(nk == NULL) continue;
				pair = next_ref_pairv(nk->pairs);
				pair->seqoff  = offset;
				pair->seqlen1 = get_tabs_len(tabs1, 9);
				pair->seqlen2 = get_tabs_len(tabs2, 9);
				pair->chr     = chr1;
				pair->lr      = rmask ^ ((flag2 & 0x10)? 1:0);
				pair->dir1    = rmask;
				pair->dir2    = 1 ^ pair->lr ^ ((flag1 & 0x10)? 1:0);
				pair->dir     = pair->dir2;
				pair->novo    = 1;
				pair->pos     = atol(get_tabs_str(tabs1, 3));
				pair->off     = rmask? len - i - 1 : i + 1 - ksize;
				used          = 1;
			}
		}
		if(used){
			ret ++;
			encap_u64list(bits, (offset + 1024) / 32);
			seq2bits(ref_u64list(bits, 0), offset, get_tabs_str(tabs1, 9), get_tabs_len(tabs1, 9));
			seq2bits(ref_u64list(bits, 0), offset + get_tabs_len(tabs1, 9), get_tabs_str(tabs2, 9), get_tabs_len(tabs2, 9));
			offset += get_tabs_len(tabs1, 9) + get_tabs_len(tabs2, 9);
		}
	}
	fprintf(stdout, "[%s] parsed %10u line\n", __FUNCTION__, rid);
	fflush(stdout);
	free_vec(tabs1);
	free_vec(tabs2);
	free_string(line1);
	free_string(line2);
	if(n_ins < 100) *avg_ins = 1000;
	else *avg_ins = t_ins / n_ins;
	return ret;
}

int scan_range(Range *r, RefSeq *rs, pairv *pairs, uint32_t lr, int off, uint32_t avg_ins){
	PairHit *pair;
	uint32_t i, j, n, x, y, a, b, chr, ins_var;
	ins_var = avg_ins / 2;
	x = y = n = 0;
	j = (ref_pairv(pairs, off)->chr == 0xFFFFU)? off + 1 : off;
	for(i=off;i<count_pairv(pairs);i++){
		pair = ref_pairv(pairs, i);
		if(pair->lr != lr) break;
		if(i <= j) continue;
		if(pair->chr == 0xFFFFU || pair->chr != ref_pairv(pairs, j)->chr || abs_sub(pair->pos, ref_pairv(pairs, j)->pos) > ins_var){
			if(i - j > n){ n = i - j; x = j; y = i; }
			j = (pair->chr == 0xFFFFU)? i + 1 : i;
		}
	}
	if(i - j > n){ n = i - j; x = j; y = i; }
	if(n == 0) return 0;
	chr = ref_pairv(pairs, x)->chr;
	a = ref_pairv(pairs, x)->pos;
	b = ref_pairv(pairs, y - 1)->pos;
	if(a > avg_ins) a -= avg_ins; else a = 0;
	if(b + avg_ins < get_u32list(rs->seqlens, chr)) b += avg_ins; else b = get_u32list(rs->seqlens, chr);
	r->id  = 0;
	r->chr = chr;
	r->off = a;
	r->len = b - a + 1;
	r->dir = ref_pairv(pairs, x)->dir;
	return i;
}

uint64_t load_novo_kmer_mates(nkv *novos, RefSeq *rs, FileReader *fr, uint32_t bin_size, uint32_t avg_ins){
	MateBin *mb;
	MateHit *mate;
	NovoKmer *nk;
	Range R1, R2, *r1, *r2;
	rangev *rg1, *rg2;
	Vector *tabs1, *tabs2;
	String *line1, *line2;
	VirtualString vl1, vl2;
	cuhash_t CU, *cu;
	uint64_t ret;
	uint32_t i, j, id, len, ins_var;
	uint32_t x, y, a, b, chr1, chr2, pos1, pos2, dir1, dir2;
	uint32_t rid, flag1, flag2, uniq1, uniq2;
	int m, n, used;
	mb = malloc(sizeof(MateBin));
	mb->bin_size = bin_size;
	mb->bins = malloc(sizeof(BitVec*) * count_strv(rs->names));
	mb->ranges = malloc(sizeof(rangev**) * count_strv(rs->names));
	for(i=0;i<count_strv(rs->names);i++){
		len = get_u32list(rs->seqlens, i);
		n = (len + bin_size - 1) / bin_size;
		mb->bins[i] = init_bitvec(n);
		mb->ranges[i] = NULL;
	}
	ins_var = avg_ins / 2;
	r1 = &R1; r2 = &R2;
	for(id=0;id<count_nkv(novos);id++){
		nk = ref_nkv(novos, id);
		sort_pairhits(ref_pairv(nk->pairs, 0), count_pairv(nk->pairs), NULL);
		if((n = scan_range(r1, rs, nk->pairs, 0, 0, avg_ins)) == 0 || scan_range(r2, rs, nk->pairs, 1, n, avg_ins) == 0) continue;
		r1->id = id * 2 + 0;
		r2->id = id * 2 + 1;
		x = r1->off / bin_size;
		y = (r1->off + r1->len - 1 + bin_size - 1) / bin_size;
		for(i=x;i<=y;i++) one_bitvec(mb->bins[r1->chr], i);
		x = r2->off / bin_size;
		y = (r2->off + r2->len - 1 + bin_size - 1) / bin_size;
		for(i=x;i<=y;i++) one_bitvec(mb->bins[r2->chr], i);
	}
	for(i=0;i<count_strv(rs->names);i++){
		index_bitvec(mb->bins[i]);
		len = (get_u32list(rs->seqlens, i) + bin_size - 1) / bin_size;
		n = rank_bitvec(mb->bins[i], len);
		mb->ranges[i] = malloc(sizeof(rangev*) * n);
		for(j=0;j<(uint32_t)n;j++) mb->ranges[i][j] = init_rangev(2);
	}
	for(id=0;id<count_nkv(novos);id++){
		nk = ref_nkv(novos, id);
		if((n = scan_range(r1, rs, nk->pairs, 0, 0, avg_ins)) == 0 || scan_range(r2, rs, nk->pairs, 1, n, avg_ins) == 0) continue;
		r1->id = id * 2 + 0;
		r2->id = id * 2 + 1;
		r1->dir = r1->dir ^ r2->dir;
		r2->dir = r1->dir;
		x = rank_bitvec(mb->bins[r1->chr], (r1->off / bin_size) + 1) - 1;
		y = rank_bitvec(mb->bins[r1->chr], (r1->off + r1->len - 1 + bin_size - 1) / bin_size + 1) - 1;
		for(i=x;i<=y;i++) push_rangev(mb->ranges[r1->chr][i], R1);
		x = rank_bitvec(mb->bins[r2->chr], (r2->off / bin_size) + 1) - 1;
		y = rank_bitvec(mb->bins[r2->chr], (r2->off + r2->len - 1 + bin_size - 1) / bin_size + 1) - 1;
		for(i=x;i<=y;i++) push_rangev(mb->ranges[r2->chr][i], R2);
	}
	for(i=0;i<count_strv(rs->names);i++){
		len = (get_u32list(rs->seqlens, i) + bin_size - 1) / bin_size;
		n = rank_bitvec(mb->bins[i], len);
		for(j=0;j<(uint32_t)n;j++) sort_range_by_id(ref_rangev(mb->ranges[i][j], 0), count_rangev(mb->ranges[i][j]), NULL);
	}
	line1 = init_string(1024);
	line2 = init_string(1024);
	tabs1 = init_vec(sizeof(VirtualString), 32);
	tabs2 = init_vec(sizeof(VirtualString), 32);
	ret = 0;
	rid = 0;
	while(1){
		chr1 = chr2 = 0xFFFFU;
		while((n = fread_line(line1, fr)) != -1){
			if(line1->string[0] == '#') continue;
			vl1.string = line1->string;
			vl1.size   = line1->size;
			clear_vec(tabs1);
			n = split_vstring(&vl1, '\t', tabs1, 1);
			if(n < 11) continue;
			flag1 = atoi(get_tabs_str(tabs1, 1));
			uniq1 = 1;
			if(flag1 & 0x04) uniq1 = 0;
			else {
				for(m=11;m<n;m++){
					if(get_tabs_len(tabs1, m) < 6) continue;
					if(memcmp(get_tabs_str(tabs1, m), "XT:A:R", 6) == 0){ uniq1 = 2; break; }
				}
			}
			if(uniq1 == 1){
				CU.key = get_tabs_str(tabs1, 2);
				CU.val = 0;
				cu = get_cuhash(rs->index, CU);
				if(cu == NULL) chr1 = 0xFFFFU;
				else chr1 = cu->val;
			} else chr1 = 0xFFFFU;
			rid ++;
			break;
		}
		if(n == -1) break;
		while((n = fread_line(line2, fr)) != -1){
			if(line2->string[0] == '#') continue;
			vl2.string = line2->string;
			vl2.size   = line2->size;
			clear_vec(tabs2);
			n = split_vstring(&vl2, '\t', tabs2, 1);
			if(n < 11) continue;
			flag2 = atoi(get_tabs_str(tabs2, 1));
			uniq2 = 1;
			if(flag2 & 0x04) uniq2 = 0;
			else {
				for(m=11;m<n;m++){
					if(get_tabs_len(tabs2, m) < 6) continue;
					if(memcmp(get_tabs_str(tabs2, m), "XT:A:R", 6) == 0){ uniq2 = 2; break; }
				}
			}
			if(uniq2 == 1){
				CU.key = get_tabs_str(tabs2, 2);
				CU.val = 0;
				cu = get_cuhash(rs->index, CU);
				if(cu == NULL) chr2 = 0xFFFFU;
				else chr2 = cu->val;
			} else chr2 = 0xFFFFU;
			rid ++;
			break;
		}
		if(n == -1) break;
		if(chr1 == 0xFFFFU || chr2 == 0xFFFFU) continue;
		pos1 = atol(get_tabs_str(tabs1, 3));
		pos2 = atol(get_tabs_str(tabs2, 3));
		dir1 = flag1 & 0x10;
		dir2 = flag2 & 0x10;
		a = pos1 / bin_size;
		b = pos2 / bin_size;
		if(get_bitvec(mb->bins[chr1], a) == 0) continue;
		if(get_bitvec(mb->bins[chr2], b) == 0) continue;
		rg1 = mb->ranges[chr1][rank_bitvec(mb->bins[chr1], a + 1) - 1];
		rg2 = mb->ranges[chr2][rank_bitvec(mb->bins[chr2], b + 1) - 1];
		a = 0; b = 0;
		used = 0;
		while(a < count_rangev(rg1) && b < count_rangev(rg2)){
			r1 = ref_rangev(rg1, a);
			r2 = ref_rangev(rg2, b);
			if((r1->id >> 1) == (r2->id >> 1)){
				a ++; b ++;
				if(r1->dir != (dir1 ^ dir2)) continue;
				nk = ref_nkv(novos, r1->id >> 1);
				mate = next_ref_matev(nk->mates);
				mate->chr1 = chr1;
				mate->chr2 = chr2;
				mate->pos1 = pos1;
				mate->pos2 = pos2;
				mate->dir1 = dir1;
				mate->dir2 = dir2;
				used = 1;
			} else if(r1->id < r2->id) a ++;
			else b ++;
		}
		if(used) ret ++;
	}
	free_string(line1);
	free_string(line2);
	free_vec(tabs1);
	free_vec(tabs2);
	return ret;
}

uint64_t output_novo_kmers(novohash *hash, RefSeq *rs, FileReader *fr, uint32_t ksize, FILE *out){
	NovoKmer *nk;
	nkv *novos;
	PairHit *pair;
	MateHit *mate;
	u64list *bits;
	uint64_t ret;
	uint32_t id, i, j, len;
	uint32_t max_l, max_r, avg_ins, bin_size;
	char *chr, seqs1[256], seqs2[256], chs[32];
	bits = init_u64list(1024);
	reset_filereader(fr);
	ret = load_novo_kmer_pairs(hash, rs, fr, ksize, bits, &avg_ins);
	fprintf(stdout, "%llu pairs contain novokmer\n", (unsigned long long)ret);
	fflush(stdout);
	bin_size = avg_ins;
	reset_filereader(fr);
	novos = init_nkv(count_novohash(hash));
	reset_iter_novohash(hash);
	while((nk = ref_iter_novohash(hash))){ push_nkv(novos, *nk); }
	free_novohash(hash);
	ret = load_novo_kmer_mates(novos, rs, fr, bin_size, avg_ins);
	fprintf(stdout, "%llu pairs locate near break points\n", (unsigned long long)ret);
	fflush(stdout);
	for(id=0;id<count_nkv(novos);id++){
		nk = ref_nkv(novos, id);
		for(i=0;i<ksize;i++) chs[ksize-1-i] = bit_base_table[(nk->kmer >> (i << 1)) & 0x03];
		chs[i] = 0;
		fprintf(out, ">%s\t%u\t%u\n", chs, count_pairv(nk->pairs), count_matev(nk->mates));
		max_l = max_r = 0;
		for(i=0;i<count_pairv(nk->pairs);i++){
			pair = ref_pairv(nk->pairs, i);
			len = pair->novo? pair->seqlen2 : pair->seqlen1;
			if(pair->off > max_l) max_l = pair->off;
			if(len - pair->off > max_r) max_r = len - pair->off;
		}
		for(i=0;i<count_pairv(nk->pairs);i++){
			pair = ref_pairv(nk->pairs, i);
			if(pair->chr != 0xFFFFU) chr = get_strv(rs->names, pair->chr)->string;
			else chr = "*";
			if(pair->novo){
				if(pair->dir1) bits2revseq(seqs1, ref_u64list(bits, 0), pair->seqoff + pair->seqlen1, pair->seqlen2);
				else bits2seq(seqs1, ref_u64list(bits, 0), pair->seqoff + pair->seqlen1, pair->seqlen2);
				if(pair->dir2) bits2revseq(seqs2, ref_u64list(bits, 0), pair->seqoff, pair->seqlen1);
				else bits2seq(seqs2, ref_u64list(bits, 0), pair->seqoff, pair->seqlen1);
			} else {
				if(pair->dir1) bits2revseq(seqs1, ref_u64list(bits, 0), pair->seqoff, pair->seqlen1);
				else bits2seq(seqs1, ref_u64list(bits, 0), pair->seqoff, pair->seqlen1);
				if(pair->dir2) bits2revseq(seqs2, ref_u64list(bits, 0), pair->seqoff + pair->seqlen1, pair->seqlen2);
				else bits2seq(seqs2, ref_u64list(bits, 0), pair->seqoff + pair->seqlen1, pair->seqlen2);
			}
			//L/R\tchr\tstrand\tpos\tnext_seq\tfirst_seq\n
			//fprintf(out, "%c\t%s\t%c\t%u\t%s\t%s\n", "RL"[pair->lr], chr, "+-"[pair->dir], pair->pos, seqs2, seqs1);
			fprintf(out, "%c\t%s\t%c\t%u\t%s\t", "RL"[pair->lr], chr, "+-"[pair->dir], pair->pos, seqs2);
			for(j=pair->off;j<max_l;j++) fprintf(out, "%c", '-');
			fprintf(out, "%s", seqs1);
			len = pair->novo? pair->seqlen2 : pair->seqlen1;
			for(j=len-pair->off;j<max_r;j++) fprintf(out, "%c", '-');
			fprintf(out, "\n");
		}
		for(i=0;i<count_matev(nk->mates);i++){
			mate = ref_matev(nk->mates, i);
			fprintf(out, "M\t%s\t%c\t%u\t%s\t%c\t%u\n", get_strv(rs->names, mate->chr1)->string, "+-"[mate->dir1], mate->pos1, 
				get_strv(rs->names, mate->chr2)->string, "+-"[mate->dir2], mate->pos2);
		}
		free_pairv(nk->pairs);
		free_matev(nk->mates);
	}
	free_u64list(bits);
	free_nkv(novos);
	return ret;
}

int usage(){
	printf("novokmer -- detecting breakpoints of SV by novo kmers in resequencing\n"
			"Author : Jue Ruan <ruanjue@gmail.com>\n"
			"Version: 1.01 (r20110926)\n"
			"Usage: novokmer -r <ref_seq_fasta_file> -i tbd.sam -o tbd.novokmer [options]\n"
			"Options:\n"
			" -r <string> Reference sequences file in fasta format *\n"
			" -m <string> Known SNP file, containing a list of records\n"
			"             \"<chr>\\t<pos>\\t<mut_base1>...\"\n"
			" -i <string> SAM file to be detected. This file should be sorted by read pairs *\n"
			" -c <string> Control SAM file.\n"
			" -o <string> Output file\n"
			" -k <int>    Kmer size, <= 27 [25]\n"
			" -a <int>    Min frequency of novo kmers [4]\n"
			" -b <int>    Max frequency of novo kmers [1000]\n"
			"\n"
			"Example1:\n"
			" novokmer -k 27 -r hg18.fa -m snps.list -i cancer.sam -c normal.sam -o cancer.novokmer -a 4\n"
			"\n"
			);
	return 1;
}

int main(int argc, char **argv){
	kmerhash *khash;
	novohash *nhash;
	RefSeq *rs;
	FileReader *inf, *ref, *mut, *ctl;
	FILE *out;
	char *infile, *outfile, *mutfile, *reffile, *ctlfile;
	int c, a, b, k;
	uint64_t ret;
	k = 25;
	a = 4;
	b = 1000;
	infile = outfile = mutfile = reffile = ctlfile = NULL;
	while((c = getopt(argc, argv, "hr:m:i:c:o:k:a:b:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'r': reffile = optarg; break;
			case 'm': mutfile = optarg; break;
			case 'i': infile = optarg; break;
			case 'c': ctlfile = optarg; break;
			case 'o': outfile = optarg; break;
			case 'k': k = atoi(optarg); break;
			case 'a': a = atoi(optarg); break;
			case 'b': b = atoi(optarg); break;
			default: return usage();
		}
	}
	if(reffile == NULL || infile == NULL || outfile == NULL) return usage();
	if(k > 27) return usage();
	if(a > 1000) return usage();
	if((inf = fopen_filereader(infile)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", infile, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	if((ref = fopen_filereader(reffile)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", reffile, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	if((out = fopen(outfile, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", outfile, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	if(mutfile){
		if((mut = fopen_filereader(mutfile)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", mutfile, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			abort();
		}
	} else mut = NULL;
	if(ctlfile){
		if((ctl = fopen_filereader(ctlfile)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", ctlfile, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			abort();
		}
	} else ctl = NULL;
	khash = build_kmerhash(inf, k);
	fprintf(stdout, "Load %llu kmer from %s\n", (unsigned long long)count_kmerhash(khash), infile);
	fflush(stdout);
	rs = build_refseq(ref, mut, k);
	fclose_filereader(ref);
	if(mut) fclose_filereader(mut);
	ret = filter_ref_genome(khash, rs, k);
	fprintf(stdout, "Filtered %llu kmer from reference\n", (unsigned long long)ret);
	fflush(stdout);
	ret = filter_separate_dbSNP(khash, rs, k);
	fprintf(stdout, "Filtered %llu kmer from separate SNP\n", (unsigned long long)ret);
	fflush(stdout);
	ret = filter_combine_dbSNP(khash, rs, k);
	fprintf(stdout, "Filtered %llu kmer from combined SNP\n", (unsigned long long)ret);
	fflush(stdout);
	if(ctl){
		ret = filter_control_kmers(khash, ctl, k);
		fprintf(stdout, "Filtered %llu kmer from Control SAM file\n", (unsigned long long)ret);
		fflush(stdout);
	}
	fprintf(stdout, "There are still %llu kmers left\n", (unsigned long long)count_kmerhash(khash));
	fflush(stdout);
	nhash = kmerhash2novohash(khash, a, b);
	fprintf(stdout, "There are %llu novokmers (occ >= %u and <= %u)\n", (unsigned long long)count_novohash(nhash), a, b);
	fflush(stdout);
	reset_filereader(inf);
	ret = output_novo_kmers(nhash, rs, inf, k, out);
	fclose(out);
	free_refseq(rs);
	fprintf(stdout, "Finished\n");
	fflush(stdout);
	return 0;
}

