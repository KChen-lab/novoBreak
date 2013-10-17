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
#include "simp_asm.h"
#include "timer.h"
#include "samtools/sam.h"
#include <stdint.h>
#include <unistd.h>

#define UNIQ_KMER_MAX_KSIZE 27
#define UNIQ_KMER_MAX_CNT 1023

#define MAX_NM 1

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
	uint64_t novo:1, lr:1, pos:32, off:14, closed:1;
} PairHit;

define_list_core(pairv, PairHit, uint32_t, 0xFU);

typedef struct {
	uint32_t chr1:16, chr2:16;
	uint32_t pos1:31, dir1:1, pos2:31, dir2:1;
} MateHit;

define_list_core(matev, MateHit, uint32_t, 0xFU);
define_list(matevec, MateHit);

typedef struct {
	uint32_t x, y;
	uint32_t lr, chr, pos, len;
} Block;

define_list(blockv, Block);

typedef struct {
	u32list *ids;
	uint32_t chr1:16, chr2:16;
	uint32_t x1, y1, x2, y2;
} BreakArea;

define_list(bkav, BreakArea);

typedef struct {
	uint64_t kmer:54, cnt:10;
	uint64_t chr:32, pos:32;
	pairv *pairs;
	matev *mates;
} NovoKmer;

define_list(nkv, NovoKmer);

#define novo_hashcode(k) u64hashcode((k).kmer)
#define novo_equals(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(novohash, NovoKmer, novo_hashcode, novo_equals);

typedef struct {
	uint32_t chr1:16, chr2:16;
	uint32_t x1, y1, x2, y2;
	pairv *pairs;
	matev *mates;
} Break;

define_list(brkv, Break);

int check_bam_alignment(bam1_t *b){
	uint32_t *cigars;
	uint8_t  *nm_aux;
	uint32_t i, n_cigar, nm;
	cigars = bam1_cigar(b);
	n_cigar = b->core.n_cigar;
	for(i=0;i<n_cigar;i++){ if((cigars[i] & 0xFU) != 0) return 0; }
	nm_aux = bam_aux_get(b, "NM");
	if(nm_aux && nm_aux[0] == 'i'){
		nm = bam_aux2i(nm_aux);
		if(nm <= MAX_NM) return 1;
	}
	return 0;
}

int check_bam_header(samfile_t *inf, RefSeq *rs){
	bam_header_t *header;
	int i;
	header = inf->header;
	if(header->n_targets != (int)count_strv(rs->names)) return 0;
	for(i=0;i<header->n_targets;i++){
		if(strcmp(header->target_name[i], get_strv(rs->names, i)->string) != 0) return 0;
		if(header->target_len[i] != get_u32list(rs->seqlens, i)) return 0;
	}
	return 1;
}

kmerhash* build_kmerhash(samfile_t *inf, uint32_t ksize){
	bam1_t *b, B;
	kmerhash *hash;
	Kmer KMER, *kmer;
	uint64_t k, r, kmask;
	uint32_t i, len, f, flag, rid;
	uint8_t *seq;
	int exists;
	hash = init_kmerhash(1023);
	KMER.cnt = 0;
	KMER.kmer = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	flag = 0;
	rid = 0;
	b = &B;
	memset(b, 0, sizeof(bam1_t));
	while(samread(inf, b) > 0){
		f = b->core.flag;
		if(flag == 0){
			if((f & 0x40) == 0) flag = 3;
			else flag = 1;
		} else if(flag == 1){
			if((f & 0x80) == 0) flag = 3;
			else flag = 0;
		}
		if(flag == 3){
			fprintf(stderr, " Two records in BAM file are not properly paired, please check your BAM file -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
			fprintf(stderr, "\tAt %s (record %u)\n", bam1_qname(b), rid + 1);
			fflush(stderr);
			abort();
		}
		if((rid & 0xFFFFU) == 0){
			fprintf(stdout, "[%s] parsed %10u alignments\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		rid ++;
		if(check_bam_alignment(b)) continue;
		//for(i=0;i<(uint32_t)b->core.l_qseq;i++) fprintf(stderr, "%c", bit4_base_table[bam1_seqi(bam1_seq(b), i)]);
		//fprintf(stderr, "\n");
		seq = bam1_seq(b);
		len = b->core.l_qseq;
		k = 0;
		for(i=0;i<len;i++){
			k = ((k << 2) | bit4_bit_table[bam1_seqi(seq, i)]) & kmask;
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
	fprintf(stdout, "[%s] parsed %10u alignments\n", __FUNCTION__, rid);
	fflush(stdout);
	if(b->data) free(b->data);
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

uint64_t filter_control_kmers(kmerhash *hash, samfile_t *ctl, uint32_t ksize){
	bam1_t *b, B;
	Kmer KMER;
	uint8_t *seq;
	uint64_t k, r, kmask, ret;
	uint32_t i, len, rid;
	KMER.cnt = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	ret = 0;
	rid = 0;
	b = &B;
	while(samread(ctl, b) > 0){
		if((rid & 0xFFFFU) == 0){
			fprintf(stdout, "[%s] parsed %10u line\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		rid ++;
		if(check_bam_alignment(b)) continue;
		seq = bam1_seq(b);
		len = b->core.l_qseq;
		k = 0;
		for(i=0;i<len;i++){
			k = ((k << 2) | bit4_bit_table[bam1_seqi(seq, i)]) & kmask;
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
	NK.chr = 0xFFFFU;
	NK.pos = 0;
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
	cmp_2nums_proc(h1.lr, h2.lr);
	cmp_2nums_proc(h1.chr, h2.chr);
	cmp_2nums_proc(h1.pos, h2.pos);
	cmp_2nums_proc(h1.seqoff, h2.seqoff);
	return 0;
	obj = obj;
}

define_quick_sort(sort_pairhits, PairHit, cmp_pairhit);

#define cmp_matehit_p1_func(h1, h2, obj) (int)(((((int64_t)(h1).chr1) << 31) | (h1).pos1) - (((int64_t)(h2).chr1 << 31) | (h2).pos1))
define_quick_sort(sort_matehits_p1, MateHit, cmp_matehit_p1_func);
define_search_array(search_matehits_p1, MateHit, cmp_matehit_p1_func);

#define cmp_matehit_p2_func(h1, h2, obj) (int)(((((int64_t)((MateHit*)obj)[h1].chr2) << 31) | ((MateHit*)obj)[h1].pos2) - ((((int64_t)((MateHit*)obj)[h2].chr2) << 31) | ((MateHit*)obj)[h2].pos2))
define_quick_sort(sort_matehits_p2, uint32_t, cmp_matehit_p2_func);
define_search_array(search_matehits_p2, uint32_t, cmp_matehit_p2_func);

int cmp_matehit_p3_func(MateHit h1, MateHit h2, void *obj){
	cmp_2nums_proc(h1.chr1, h2.chr1);
	cmp_2nums_proc(h1.chr2, h2.chr2);
	cmp_2nums_proc(h1.pos1, h2.pos1);
	cmp_2nums_proc(h1.pos2, h2.pos2);
	return 0;
	obj = obj;
}

define_quick_sort(sort_matehits, MateHit, cmp_matehit_p3_func);

#define cmp_nk_breakpoint(n1, n2, obj) ((((int64_t)(((int64_t)(n1).chr) << 32) | (n1).pos)) - ((int64_t)((((int64_t)(n2).chr) << 32) | (n2).pos)))
define_quick_sort(sort_novokmers_by_bp, NovoKmer, cmp_nk_breakpoint);

int cmp_bka_func(BreakArea b1, BreakArea b2, void *obj){
	cmp_2nums_proc(b1.chr1, b2.chr1);
	cmp_2nums_proc(b1.x1, b2.x1);
	return 0;
	obj = obj;
}

define_quick_sort(sort_bkas, BreakArea, cmp_bka_func);

#define get_tabs_str(tabs, idx) ((VirtualString*)get_vec_ref(tabs, idx))->string
#define get_tabs_len(tabs, idx) ((VirtualString*)get_vec_ref(tabs, idx))->size

inline uint32_t abs_sub(uint32_t a, uint32_t b){ return (a < b)? b - a : a - b; }

int read_bam_alignment(samfile_t *fr, bam1_t *b, uint32_t *flag, uint32_t *chr, uint32_t *pos, char *seq, uint32_t *len){
	uint8_t *aux;
	int uniq, i;
	while(samread(fr, b) > 0){
		*flag = b->core.flag;
		uniq = 1;
		if((*flag) & 0x04) uniq = 0;
		else if((aux = bam_aux_get(b, "XT"))){
			if(bam_aux2A(aux) == 'R') uniq = 0;
		}
		if(uniq == 1){
			if(b->core.tid < 0) *chr = 0xFFFFU;
			else *chr = b->core.tid;
		} else *chr = 0xFFFFU;
		*pos = b->core.pos;
		if(len) *len = b->core.l_qseq;
		if(seq){
			aux = bam1_seq(b);
			for(i=0;i<b->core.l_qseq;i++) seq[i] = bit4_base_table[bam1_seqi(aux, i)];
			seq[i] = '\0';
		}
		return 1;
	}
	return 0;
}

uint64_t load_novo_kmer_pairs(novohash *hash, samfile_t *inf, uint32_t ksize, u64list *bits, uint32_t *avg_ins){
	bam1_t *b1, *b2, B1, B2;
	NovoKmer NK, *nk;
	PairHit *pair;
	uint64_t ret, offset, k, r, kmask, t_ins, n_ins;
	uint32_t i, n, len, rid, rmask, used, flag1, flag2, chr1, chr2, pos1, pos2, len1, len2;
	char *seq, seq1[256], seq2[256];
	ret = 0;
	NK.pairs = NULL;
	NK.mates = NULL;
	NK.cnt = 0;
	NK.chr = 0xFFFFU;
	NK.pos = 0;
	kmask = (1LLU << (2 * ksize)) - 1;
	offset = 0;
	rid = 0;
	t_ins = 0;
	n_ins = 0;
	b1 = &B1; b2 = &B2;
	memset(b1, 0, sizeof(bam1_t));
	memset(b2, 0, sizeof(bam1_t));
	while(1){
		if(read_bam_alignment(inf, b1, &flag1, &chr1, &pos1, seq1, &len1) == 0) break;
		if(read_bam_alignment(inf, b2, &flag2, &chr2, &pos2, seq2, &len2) == 0) break;
		if((rid & 0xFFFFU) == 0){
			fprintf(stdout, "[%s] parsed %10u pairs\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		rid ++;
		n = abs_sub(pos1, pos2);
		if(chr1 == chr2 && chr1 != 0xFFFFU && n < 10 * 1024){
			t_ins += n + len1;
			n_ins ++;
		}
		used = 0;
		for(n=0;n<1;n++){
			if(check_bam_alignment(b1)) continue;
			seq = seq1; len = len1;
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
//				fprintf(stdout, "\nchr1:%u\tflag:%u\trmask:%u\tpos1:%u\t%s\n", chr1, flag1, rmask, pos1, seq);
//				fflush(stdout);
				pair = next_ref_pairv(nk->pairs);
				pair->closed  = 1;
				pair->seqoff  = offset;
				pair->seqlen1 = len1;
				pair->seqlen2 = len2;
				pair->chr     = chr2;
				pair->lr      = rmask ^ ((flag1 & 0x10)? 1:0); 
				pair->dir1    = rmask;
				pair->dir2    = 1 ^ pair->lr ^ ((flag2 & 0x10)? 1:0);
				pair->dir     = pair->dir2;
				pair->novo    = 0;
				pair->pos     = pos2;
				pair->off     = rmask? len - i - 1 : i + 1 - ksize;
				used          = 1;
//				fprintf(stdout, "chr:%u\tlr:%u\tdir1:%u\tdir2:%u\tdir:%u\tpos:%u\n", pair->chr, pair->lr, pair->dir1, pair->dir2, pair->dir, pair->pos);
//				fflush(stdout);
			}
		}
		for(n=0;n<1;n++){
			if(check_bam_alignment(b2)) continue;
			seq = seq2; len = len2;
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
//				fprintf(stdout, "\nchr2:%u\tflag:%u\trmask:%u\tpos2:%u\t%s\n", chr2, flag2, rmask, pos2, seq);
//				fflush(stdout);
				pair = next_ref_pairv(nk->pairs);
				pair->closed  = 1;
				pair->seqoff  = offset;
				pair->seqlen1 = len1;
				pair->seqlen2 = len2;
				pair->chr     = chr1;
				pair->lr      = rmask ^ ((flag2 & 0x10)? 1:0);
				pair->dir1    = rmask;
				pair->dir2    = 1 ^ pair->lr ^ ((flag1 & 0x10)? 1:0);
				pair->dir     = pair->dir2;
				pair->novo    = 1;
				pair->pos     = pos1;
				pair->off     = rmask? len - i - 1 : i + 1 - ksize;
				used          = 1;
//				fprintf(stdout, "chr:%u\tlr:%u\tdir1:%u\tdir2:%u\tdir:%u\tpos:%u\n", pair->chr, pair->lr, pair->dir1, pair->dir2, pair->dir, pair->pos);
//				fflush(stdout);
			}
		}
		if(used){
			ret ++;
			encap_u64list(bits, (offset + 1024) / 32);
			seq2bits(ref_u64list(bits, 0), offset, seq1, len1);
			seq2bits(ref_u64list(bits, 0), offset + len1, seq2, len2);
			offset += len1 + len2;
		}
	}
	fprintf(stdout, "[%s] parsed %10u pairs\n", __FUNCTION__, rid);
	fflush(stdout);
	if(b1->data) free(b1->data);
	if(b2->data) free(b2->data);
	if(n_ins < 100) *avg_ins = 1000;
	else *avg_ins = t_ins / n_ins;
	return ret;
}

uint32_t search_blocks(matevec *mates, u32list *mids, Block *b1, Block *b2, matev *rs){
	MateHit MATE, *mate;
	Block *t;
	int64_t idx;
	uint32_t k, ranges[4], cmp, ret;
	if(b1->chr < b2->chr) cmp = 0;
	else if(b1->chr > b2->chr) cmp = 1;
	else if(b1->pos < b2->pos) cmp = 0;
	else cmp = 1;
	if(cmp){ swap_tmp(b1, b2, t); }
	MATE.dir1 = MATE.dir2 = 0;
	MATE.chr1 = b1->chr;
	MATE.pos1 = b1->pos;
	idx = search_matehits_p1(ref_matevec(mates, 0), count_matevec(mates), MATE, NULL);
	if(idx < 0) idx = - 1 - idx;
	ranges[0] = idx;
	MATE.chr1 = b1->chr;
	MATE.pos1 = b1->pos + b1->len - 1;
	idx = search_matehits_p1(ref_matevec(mates, 0), count_matevec(mates), MATE, NULL);
	if(idx < 0) idx = - 1 - idx;
	else idx ++;
	ranges[1] = idx;
	if(ranges[1] <= ranges[0]) return 0;
	MATE.chr2 = b2->chr;
	MATE.pos2 = b2->pos;
	push_u32list(mids, count_matevec(mates));
	trunc_u32list(mids, 1);
	push_matevec(mates, MATE);
	trunc_matevec(mates, 1);
	idx = search_matehits_p2(ref_u32list(mids, 0), count_u32list(mids), count_matevec(mates), ref_matevec(mates, 0));
	if(idx < 0) idx = - 1 - idx;
	ranges[2] = idx;
	MATE.chr2 = b2->chr;
	MATE.pos2 = b2->pos + b2->len - 1;
	push_u32list(mids, count_matevec(mates));
	trunc_u32list(mids, 1);
	push_matevec(mates, MATE);
	trunc_matevec(mates, 1);
	idx = search_matehits_p2(ref_u32list(mids, 0), count_u32list(mids), count_matevec(mates), ref_matevec(mates, 0));
	if(idx < 0) idx = - 1 - idx;
	else idx ++;
	ranges[3] = idx;
	if(ranges[3] <= ranges[2]) return 0;
	ret = 0;
	for(k=ranges[2];k<ranges[3];k++){
		if(get_u32list(mids, k) < ranges[1] && get_u32list(mids, k) >= ranges[0]){
			mate = ref_matevec(mates, get_u32list(mids, k));
			if(mate->chr1 != b1->chr) continue;
			if(mate->pos1 < b1->pos) continue;
			if(mate->pos1 >= b1->pos + b1->len) continue;
			if(mate->chr2 != b2->chr) continue;
			if(mate->pos2 < b2->pos) continue;
			if(mate->pos2 >= b2->pos + b2->len) continue;
			//fprintf(stdout, "M\t%s\t%c\t%u\t%s\t%c\t%u\n", get_strv(rs->names, mate->chr1)->string, "+-"[mate->dir1], mate->pos1, get_strv(rs->names, mate->chr2)->string, "+-"[mate->dir2], mate->pos2);
			push_matev(rs, *mate);
			ret ++;
		}
	}
	return ret;
}

uint64_t load_novo_kmer_mates(nkv *novos, RefSeq *rs, samfile_t *inf, uint32_t avg_ins){
	bam1_t *b1, *b2, B1, B2;
	PairHit *pair;
	MateHit *mate;
	matevec *mates;
	u32list *mids;
	NovoKmer *nk;
	blockv *bs;
	Block *block;
	uint64_t ret;
	uint32_t i, j, id, ins_var, mid;
	uint32_t x, y, chr1, chr2, pos1, pos2, dir1, dir2, tmp;
	uint32_t rid, flag1, flag2;
	ins_var = avg_ins;
	mates = init_matevec(1024);
	ret = 0;
	rid = 0;
	b1 = &B1; b2 = &B2;
	memset(b1, 0, sizeof(bam1_t));
	memset(b2, 0, sizeof(bam1_t));
	while(1){
		if(read_bam_alignment(inf, b1, &flag1, &chr1, &pos1, NULL, NULL) == 0) break;
		if(read_bam_alignment(inf, b2, &flag2, &chr2, &pos2, NULL, NULL) == 0) break;
		if((rid & 0xFFFFU) == 0){
			fprintf(stdout, "[%s] parsed %10u pairs\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		rid ++;
		if(chr1 == 0xFFFFU || chr2 == 0xFFFFU) continue;
		dir1 = (flag1 & 0x10)? 1 : 0;
		dir2 = (flag2 & 0x10)? 1 : 0;
		if(dir1 != dir2 && chr1 == chr2 && abs_sub(pos1, pos2) < 2 * avg_ins) continue;
		//fprintf(stderr, "%s\t%c\t%u\t%s\t%c\t%u\n", get_strv(rs->names, chr1)->string, "+-"[dir1], pos1, get_strv(rs->names, chr2)->string, "+-"[dir2], pos2);
		mate = next_ref_matevec(mates);
		if(chr1 > chr2 || (chr1 == chr2 && pos1 > pos2)){
			swap_tmp(chr1, chr2, tmp);
			swap_tmp(pos1, pos2, tmp);
			dir1 = !dir1; dir2 = !dir2;
			swap_tmp(dir1, dir2, tmp);
		}
		mate->chr1 = chr1;
		mate->chr2 = chr2;
		mate->dir1 = dir1;
		mate->dir2 = dir2;
		mate->pos1 = pos1;
		mate->pos2 = pos2;
	}
	fprintf(stdout, "[%s] found %10u abnormal pairs\n", __FUNCTION__, (unsigned)count_matevec(mates));
	fflush(stdout);
	if(b1->data) free(b1->data);
	if(b2->data) free(b2->data);
	sort_matehits_p1(ref_matevec(mates, 0), count_matevec(mates), NULL);
	mids  = init_u32list(count_matevec(mates));
	for(i=0;i<count_matevec(mates);i++) push_u32list(mids, i);
	sort_matehits_p2(ref_u32list(mids, 0), count_u32list(mids), ref_matevec(mates, 0));
	fprintf(stdout, "[%s] indexed abnormal pairs\n", __FUNCTION__);
	fflush(stdout);
	bs = init_blockv(32);
	ret = 0;
	for(id=0;id<count_nkv(novos);id++){
		nk = ref_nkv(novos, id);
		mid = 0;
		clear_blockv(bs);
		block = next_ref_blockv(bs);
		block->x   = 0;
		block->lr  = 0;
		block->chr = 0xFFFFU;
		block->pos = 0;
		//fprintf(stdout, ">%u\n", id);
		for(i=0;i<count_pairv(nk->pairs);i++){
			pair = ref_pairv(nk->pairs, i);
			//fprintf(stdout, "%c\t%u\t%c\t%u\n", "LR"[pair->lr], pair->chr, "+-"[pair->dir], (unsigned)pair->pos);
			if(pair->lr != block->lr || pair->chr != block->chr || abs_sub(pair->pos, block->pos) > ins_var){
				if(pair->lr != block->lr){ mid = (block->chr == 0xFFFFU)? count_blockv(bs) - 1 : count_blockv(bs); }
				if(block->chr != 0xFFFFU){ block->y = i; block = next_ref_blockv(bs); }
				block->lr  = pair->lr;
				block->chr = pair->chr;
				block->x   = i;
				block->pos = pair->pos;
			}
		}
		//fprintf(stdout, "\n");
		if(block->chr != 0xFFFFU) block->y = i;
		else trunc_blockv(bs, 1);
		if(mid == 0 || mid  == count_blockv(bs)) continue;
		nk->chr = ref_blockv(bs, 0)->chr;
		nk->pos = ref_blockv(bs, 0)->pos;
		if(ref_blockv(bs, mid)->chr < nk->chr || (ref_blockv(bs, mid)->chr == nk->chr && ref_blockv(bs, mid)->pos < nk->pos)){
			nk->chr = ref_blockv(bs, mid)->chr;
			nk->pos = ref_blockv(bs, mid)->pos;
		}
		for(i=0;i<count_blockv(bs);i++){
			block = ref_blockv(bs, i);
			x = ref_pairv(nk->pairs, block->x)->pos;
			y = ref_pairv(nk->pairs, block->y - 1)->pos;
			if(x > 2 * avg_ins) x -= 2 * avg_ins;
			else x = 0;
			if(y + 2 * avg_ins < get_u32list(rs->seqlens, block->chr)) y += 2 * avg_ins;
			else y = get_u32list(rs->seqlens, block->chr) - 1;
			block->pos = x; block->len = y - x + 1;
			//fprintf(stdout, "%c\t%u\t%u\t%u\n", "RL"[block->lr], block->chr, block->pos, block->len);
		}
		//fprintf(stdout, "\n");
		//fprintf(stdout, " -- mid = %u in %s -- %s:%d --\n", mid, __FUNCTION__, __FILE__, __LINE__);
		x = 0; y = 0;
		for(i=0;i<mid;i++){
			block = ref_blockv(bs, i);
			if(block->y - block->x > y){ y = block->y - block->x; x = i; }
		}
		i = x;
		x = mid; y = 0;
		for(j=mid;j<count_blockv(bs);j++){
			block = ref_blockv(bs, j);
			if(block->y - block->x > y){ y = block->y - block->x; x = j; }
		}
		j = x;
		{
			block = ref_blockv(bs, i);
			for(x=block->x;x<block->y;x++) ref_pairv(nk->pairs, x)->closed = 0;
			block = ref_blockv(bs, j);
			for(x=block->x;x<block->y;x++) ref_pairv(nk->pairs, x)->closed = 0;
		}
		if(ref_blockv(bs, i)->chr != ref_blockv(bs, j)->chr || (
					ref_blockv(bs, i)->pos + ref_blockv(bs, i)->len + 3 * avg_ins < ref_blockv(bs, j)->pos ||
					ref_blockv(bs, j)->pos + ref_blockv(bs, j)->len + 3 * avg_ins < ref_blockv(bs, i)->pos)){
			ret += search_blocks(mates, mids, ref_blockv(bs, i), ref_blockv(bs, j), nk->mates);
		}
	}
	free_blockv(bs);
	free_matevec(mates);
	free_u32list(mids);
	return ret;
}

bkav* mates2bkas(nkv *novos, uint32_t ins_var){
	bkav *bkas;
	BreakArea *bka, *bka2;
	MateHit *mate;
	NovoKmer *nk;
	uint32_t id, i, f;
	bkas = init_bkav(count_nkv(novos));
	for(id=0;id<count_nkv(novos);id++){
		nk = ref_nkv(novos, id);
		bka = next_ref_bkav(bkas);
		bka->ids = init_u32list(2);
		push_u32list(bka->ids, id);
		bka->chr1 = 0xFFFFU;
		bka->chr2 = 0xFFFFU;
		bka->x1   = 0;
		bka->y1   = 0;
		bka->x2   = 0;
		bka->y2   = 0;
		if(count_matev(nk->mates) == 0) continue;
		bka->chr1 = ref_matev(nk->mates, 0)->chr1;
		bka->chr2 = ref_matev(nk->mates, 0)->chr2;
		bka->x1   = ref_matev(nk->mates, 0)->pos1;
		bka->y1   = ref_matev(nk->mates, count_matev(nk->mates) - 1)->pos1;
		bka->x2   = ref_matev(nk->mates, 0)->pos2;
		bka->y2   = ref_matev(nk->mates, count_matev(nk->mates) - 1)->pos2;
		for(i=0;i<count_matev(nk->mates);i++){
			mate = ref_matev(nk->mates, i);
			if(mate->pos1 < bka->x1) bka->x1 = mate->pos1;
			if(mate->pos1 > bka->y1) bka->y1 = mate->pos1;
		}
	}
	if(count_bkav(bkas) == 0) return bkas;
	sort_bkas(ref_bkav(bkas, 0), count_bkav(bkas), NULL);
	bka = ref_bkav(bkas, 0);
	for(i=1;i<count_bkav(bkas);i++){
		bka2 = ref_bkav(bkas, i);
		f = 0;
		if(bka->chr1 != 0xFFFFU && bka->chr2 != 0xFFFFU && bka2->chr1 == bka->chr1 && bka2->chr2 == bka->chr2){
			if(abs_sub(bka2->x1, bka->x1) <= ins_var && abs_sub(bka2->y1, bka->y1) <= ins_var){
				if(abs_sub(bka2->x2, bka->x2) <= ins_var && abs_sub(bka2->y2, bka->y2) <= ins_var) f = 1;
			}
		}
		if(f){
			append_u32list(bka->ids, bka2->ids);
			clear_u32list(bka2->ids);
			if(bka->x1 > bka2->x1) bka->x1 = bka2->x1;
			if(bka->x2 > bka2->x2) bka->x2 = bka2->x2;
			if(bka->y1 < bka2->y1) bka->y1 = bka2->y1;
			if(bka->y2 < bka2->y2) bka->y2 = bka2->y2;
		} else {
			bka = bka2;
		}
	}
	return bkas;
}

brkv* transform_novokmers(nkv *novos, bkav *bkas){
	brkv *brks;
	pairv *pairs;
	matev *mates;
	Break *brk;
	BreakArea *bka;
	NovoKmer *nk;
	PairHit *pair;
	MateHit *mate;
	uint32_t i, j, s;
	brks = init_brkv(1024);
	pairs = init_pairv(1024);
	mates = init_matev(1024);
	for(i=0;i<count_bkav(bkas);i++){
		bka = ref_bkav(bkas, i);
		if(count_u32list(bka->ids) == 0) continue;
		brk = next_ref_brkv(brks);
		brk->chr1 = bka->chr1;
		brk->chr2 = bka->chr2;
		brk->x1   = bka->x1;
		brk->y1   = bka->y1;
		brk->x2   = bka->x2;
		brk->y2   = bka->y2;
		brk->pairs = init_pairv(4);
		brk->mates = init_matev(4);
		clear_pairv(pairs);
		clear_matev(mates);
		for(j=0;j<count_u32list(bka->ids);j++){
			nk = ref_nkv(novos, get_u32list(bka->ids, j));
			for(s=0;s<count_pairv(nk->pairs);s++){ pair = ref_pairv(nk->pairs, s); if(pair->closed == 0) push_pairv(pairs, *pair); }
			append_matev(mates, nk->mates);
			free_pairv(nk->pairs);
			free_matev(nk->mates);
		}
		sort_pairhits(ref_pairv(pairs, 0), count_pairv(pairs), NULL);
		pair = NULL;
		for(j=0;j<count_pairv(pairs);j++){
			if(pair == NULL){ pair = ref_pairv(pairs, j); continue; }
			if(pair->seqoff != ref_pairv(pairs, j)->seqoff){
				push_pairv(brk->pairs, *pair);
				pair = ref_pairv(pairs, j);
			}
		}
		if(pair) push_pairv(brk->pairs, *pair);
		sort_matehits(ref_matev(mates, 0), count_matev(mates), NULL);
		mate = NULL;
		for(j=0;j<count_matev(mates);j++){
			if(mate == NULL){ mate = ref_matev(mates, j); continue; }
			if(memcmp(mate, ref_matev(mates, j), sizeof(MateHit)) != 0){
				push_matev(brk->mates, *mate);
				mate = ref_matev(mates, j);
			}
		}
		if(mate) push_matev(brk->mates, *mate);
	}
	free_pairv(pairs);
	free_matev(mates);
	free_nkv(novos);
	free_bkav(bkas);
	return brks;
}

uint64_t output_novo_kmers(RefSeq *rs, brkv *brks, uint64_t *bits, int asm_no_mate,  FILE *out){
	PairHit *pair;
	MateHit *mate;
	Break *brk;
	SimpAssembler *sa;
	SimpContigInfo *ctg;
	uint64_t ret;
	uint32_t id, i, f;
	char *chr1, *chr2, *chr, seqs1[256];
	ret = 0;
	sa = init_simpasm(5, 0, 3, 20, 0.90, 6, 1);
	for(id=0;id<count_brkv(brks);id++){
		if((id & 0xF) == 0){
			fprintf(stdout, "[%s] process %u/%u/%u breakpoints\r", __FUNCTION__, (unsigned)ret, id, (unsigned)count_brkv(brks));
			fflush(stdout);
		}
		brk = ref_brkv(brks, id);
		f = count_matev(brk->mates)? 1 : asm_no_mate;
		if(f == 0) continue;
		ret ++;
		reset_simpasm(sa);
		chr1 = (brk->chr1 == 0xFFFFU)? "*" : get_strv(rs->names, brk->chr1)->string;
		chr2 = (brk->chr2 == 0xFFFFU)? "*" : get_strv(rs->names, brk->chr2)->string;
		fprintf(out, ">break%u\t%s\t%u\t%u\t%s\t%u\t%u\t%u\t%u\n",
				id, chr1, brk->x1, brk->y1, chr2, brk->x2, brk->y2, (unsigned)count_pairv(brk->pairs), (unsigned)count_matev(brk->mates));
		for(i=0;i<count_pairv(brk->pairs);i++){
			pair = ref_pairv(brk->pairs, i);
			if(pair->chr != 0xFFFFU) chr = get_strv(rs->names, pair->chr)->string;
			else chr = "*";
			if(pair->novo){
				if(pair->dir1) bits2revseq(seqs1, bits, pair->seqoff + pair->seqlen1, pair->seqlen2);
				else bits2seq(seqs1, bits, pair->seqoff + pair->seqlen1, pair->seqlen2);
				//if(pair->dir2) bits2revseq(seqs2, ref_u64list(bits, 0), pair->seqoff, pair->seqlen1);
				//else bits2seq(seqs2, ref_u64list(bits, 0), pair->seqoff, pair->seqlen1);
			} else {
				if(pair->dir1) bits2revseq(seqs1, bits, pair->seqoff, pair->seqlen1);
				else bits2seq(seqs1, bits, pair->seqoff, pair->seqlen1);
				//if(pair->dir2) bits2revseq(seqs2, ref_u64list(bits, 0), pair->seqoff + pair->seqlen1, pair->seqlen2);
				//else bits2seq(seqs2, ref_u64list(bits, 0), pair->seqoff + pair->seqlen1, pair->seqlen2);
			}
			fprintf(out, "P\t%s\t%u\t%s\n", chr, pair->pos, seqs1);
			if(f && pair->novo) push_simpasm(sa, i, seqs1, strlen(seqs1), 1);
		}
		{
			simple_assemble(sa);
			begin_iter_simpasm(sa);
			while((ctg = iter_simpasm(sa))){
				fprintf(out, "C\t%s\n", ctg->seq->string);
			}
		}
		for(i=0;i<count_matev(brk->mates);i++){
			mate = ref_matev(brk->mates, i);
			fprintf(out, "M\t%s\t%c\t%u\t%s\t%c\t%u\n", get_strv(rs->names, mate->chr1)->string, "+-"[mate->dir1], mate->pos1, 
				get_strv(rs->names, mate->chr2)->string, "+-"[mate->dir2], mate->pos2);
		}
		free_pairv(brk->pairs);
		free_matev(brk->mates);
	}
	fprintf(stdout, "[%s] process %10u/%u breakpoints\n", __FUNCTION__, id, (unsigned)count_brkv(brks));
	fflush(stdout);
	free_simpasm(sa);
	free_brkv(brks);
	return ret;
}

int usage(){
	printf("novokmer -- detecting breakpoints of SV by novo kmers in resequencing\n"
			"Author : Jue Ruan <ruanjue@gmail.com>\n"
			"Version: 1.03 (r20111008)\n"
			"Usage: novokmer -r <ref_seq_fasta_file> -i tbd.sam -o tbd.novokmer [options]\n"
			"Options:\n"
			" -r <string> Reference sequences file in fasta format *\n"
			" -m <string> Known SNP file, containing a list of records\n"
			"             \"<chr>\\t<pos>\\t<mut_base1>...\"\n"
			" -i <string> BAM file to be detected. This file should be sorted by read pairs *\n"
			" -c <string> Control BAM file.\n"
			" -o <string> Output file\n"
			" -k <int>    Kmer size, <= 27 [25]\n"
			" -a <int>    Min frequency of novo kmers [4]\n"
			" -b <int>    Max frequency of novo kmers [1000]\n"
			" -f          Force assemble breakpoints without paired-end coverage [no]\n"
			"\n"
			"Example1:\n"
			" novokmer -k 27 -r hg18.fa -m snps.list -i cancer.bam -c normal.bam -o cancer.novokmer -a 4\n"
			"\n"
			);
	return 1;
}

int main(int argc, char **argv){
	kmerhash *khash;
	novohash *nhash;
	nkv *novos;
	bkav *bkas;
	brkv *brks;
	RefSeq *rs;
	u64list *bits;
	NovoKmer *nk;
	FileReader *ref, *mut;
	samfile_t *inf, *ctl;
	FILE *out;
	uint32_t avg_ins, ins_var;
	char *infile, *outfile, *mutfile, *reffile, *ctlfile;
	int c, a, b, k, asm_no_mate;
	uint64_t ret;
	k = 25;
	a = 4;
	b = 1000;
	asm_no_mate = 0;
	infile = outfile = mutfile = reffile = ctlfile = NULL;
	while((c = getopt(argc, argv, "hr:m:i:c:o:k:a:b:f")) != -1){
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
			case 'f': asm_no_mate = 1; break;
			default: return usage();
		}
	}
	if(reffile == NULL || infile == NULL || outfile == NULL) return usage();
	if(k > 27) return usage();
	if(a > 1000) return usage();
	if((inf = samopen(infile, "rb", NULL)) == NULL){
		fprintf(stderr, " -- Cannot open bam file %s in %s -- %s:%d --\n", infile, __FUNCTION__, __FILE__, __LINE__);
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
		if((ctl = samopen(ctlfile, "rb", NULL)) == NULL){
			fprintf(stderr, " -- Cannot open bamfile %s in %s -- %s:%d --\n", ctlfile, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			abort();
		}
	} else ctl = NULL;
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	rs = build_refseq(ref, mut, k);
	fclose_filereader(ref);
	if(mut) fclose_filereader(mut);
	fprintf(stdout, "Load reference\n");
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	if(check_bam_header(inf, rs) == 0){
		fprintf(stdout, "[Fetal Error] Reference file is inconsistent with BAM file\n");
		fflush(stdout);
		abort();
	}
	khash = build_kmerhash(inf, k);
	samclose(inf);
	fprintf(stdout, "Load %llu kmer from %s\n", (unsigned long long)count_kmerhash(khash), infile);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	ret = filter_ref_genome(khash, rs, k);
	fprintf(stdout, "Filtered %llu kmer from reference\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	ret = filter_separate_dbSNP(khash, rs, k);
	fprintf(stdout, "Filtered %llu kmer from separate SNP\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	ret = filter_combine_dbSNP(khash, rs, k);
	fprintf(stdout, "Filtered %llu kmer from combined SNP\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	if(ctl){
		ret = filter_control_kmers(khash, ctl, k);
		fprintf(stdout, "Filtered %llu kmer from Control BAM file\n", (unsigned long long)ret);
		fflush(stdout);
		fprintf(stdout, "[%s]\n", date()); fflush(stdout);
		samclose(ctl);
	}
	fprintf(stdout, "There are still %llu kmers left\n", (unsigned long long)count_kmerhash(khash));
	fflush(stdout);
	nhash = kmerhash2novohash(khash, a, b);
	fprintf(stdout, "There are %llu novokmers (occ >= %u and <= %u)\n", (unsigned long long)count_novohash(nhash), a, b);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	bits = init_u64list(1024);
	inf = samopen(infile, "rb", NULL);
	ret = load_novo_kmer_pairs(nhash, inf, k, bits, &avg_ins);
	samclose(inf);
	fprintf(stdout, "Average insert size was estimated to %u bp\n", avg_ins);
	fprintf(stdout, "%llu pairs contain novokmer\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	ins_var = avg_ins / 2;
	novos = init_nkv(count_novohash(nhash));
	reset_iter_novohash(nhash);
	while((nk = ref_iter_novohash(nhash))){
		sort_pairhits(ref_pairv(nk->pairs, 0), count_pairv(nk->pairs), NULL);
		push_nkv(novos, *nk);
	}
	free_novohash(nhash);
	inf = samopen(infile, "rb", NULL);
	ret = load_novo_kmer_mates(novos, rs, inf, avg_ins);
	samclose(inf);
	sort_novokmers_by_bp(ref_nkv(novos, 0), count_nkv(novos), NULL);
	fprintf(stdout, "%llu pairs locate near break points (may be some duplications)\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	bkas = mates2bkas(novos, ins_var);
	brks = transform_novokmers(novos, bkas);
	fprintf(stdout, "Clustered breakpoints\n");
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	ret = output_novo_kmers(rs, brks, ref_u64list(bits, 0), asm_no_mate, out);
	fprintf(stdout, "Output %llu breakpoints\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	free_u64list(bits);
	fclose(out);
	free_refseq(rs);
	fprintf(stdout, "Finished\n");
	fflush(stdout);
	return 0;
}
