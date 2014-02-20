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


#include "filter_kmer.h"

kmerhash* build_kmerhash(FileReader *fr, uint32_t ksize, int is_fq, kmerhash* hash) {
	kmer_t KMER, *kmer;
	uint64_t k, r, kmask;
	uint32_t i, j, len, rid, lowq_num;
	int exists;
	Sequence *seq;
	char *quality;
	rid = 0;
	KMER.cnt = 0;
	KMER.cnt2 = 0;
	KMER.kmer = 0;
//	kmask = (1LLU << (2 * ksize)) - 1;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	seq = NULL;
	while (is_fq?fread_fastq_adv(&seq, fr, FASTQ_FLAG_NO_NAME):fread_fasta_adv(&seq, fr, FASTA_FLAG_NO_NAME)) {
		rid ++;
		if ((rid & 0xFFFFU) == 0) {
			fprintf(stdout, "[%s] parsed %10u reads\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		k = 0;
		len = seq->seq.size;
		quality = seq->qual.string;
		lowq_num = 0;
		for (j = 0; j < len; j++) {
			if (quality[j] <= '#') lowq_num ++;
		}
		if ((float)lowq_num/len >= 0.5) continue; // too many low quality bases

		for (i = 0; i < ksize-1; i++) {
			k = (k << 2) | base_bit_table[(int)seq->seq.string[i]];
		}
		for (i = 0; i <= len-ksize; i++) {
			if (quality[i+ksize-1] <= '#' && i+ksize < len && quality[i+ksize] <= '#') continue; // trim ends with low qualities
			k = ((k << 2) | base_bit_table[(int)seq->seq.string[i+ksize-1]])  & kmask;
			//if (i + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			kmer = prepare_kmerhash(hash, KMER, &exists);
			if (exists) {
				if (kmer->cnt < UNIQ_KMER_MAX_CNT)
					kmer->cnt ++;
			} else {
				kmer->kmer = KMER.kmer;
				kmer->cnt = 1;
				kmer->cnt2 = 0;
			}
		}
	}
	fprintf(stdout, "[%s] processed %10u reads\n", __FUNCTION__, rid);
	fflush(stdout);
	return hash;
}

uint64_t filter_ref_kmers(kmerhash *hash, FileReader *fr, uint32_t ksize) {
	Sequence *seq;
	kmer_t KMER;
	uint64_t k, r, kmask, ret = 0;
	uint32_t i, len;
	seq = NULL;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	while (fread_fasta(&seq, fr)) {
		fprintf(stdout, "Filtering %s\r", seq->name.string);
		fflush(stdout);
		k = 0;
		len = seq->seq.size;
		for (i = 0; i < ksize-1; i++)
			k = (k << 2) | base_bit_table[(int)seq->seq.string[i]];
		for (i = 0; i <= len-ksize; i++) {
			k = ((k << 2) | base_bit_table[(int)seq->seq.string[i+ksize-1]])  & kmask;
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			ret += remove_kmerhash(hash, KMER);
		}
	}
	return ret;
}

pairv* loadkmerseq(kmerhash *hash, uint32_t ksize, uint32_t mincnt, uint32_t maxcnt2, FileReader *f1, FileReader *f2, int is_somatic) {
	pairv *pairs = init_pairv(4);
	pair_t *pair;
	kmer_t KMER, *ret;
	Sequence *read[2];
	int is_fq[2];
	uint32_t rid, len, i, j, ii, lowq_num; //, len_head;
	char *seq, *quality, kmer_seq[256];
	uint64_t k, kmask = (1LLU << (2 * ksize)) - 1, r;

	read[0] = read[1] = NULL;
	rid = 0;
	KMER.cnt = 0;

	is_fq[0] = (guess_seq_file_type(f1) == 2);
	is_fq[1] = (guess_seq_file_type(f2) == 2);

	while (1) {
		//label:
		if (!(is_fq[0]?fread_fastq(&read[0], f1):fread_fasta(&read[0], f1)))  break;
		if (!(is_fq[1]?fread_fastq(&read[1], f2):fread_fasta(&read[1], f2)))  break;
		rid ++;
		if ((rid & 0xFFFFU) == 0) {
			fprintf(stdout, "[%s] parsed %10u pairs\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		for (i = 0; i < 2; i ++) {
			seq = read[i]->seq.string;
			len = read[i]->seq.size;
			quality = read[i]->qual.string;
			lowq_num = 0;
			for (j = 0; j < len; j++) {
				if (quality[j] <= '#') lowq_num ++;
			}
			if ((float)lowq_num/len >= 0.5) continue; // too many low quality bases
			k = 0;
			//fprintf(stderr, "seq%d\t%s\n", i, seq);
			for (j = 0; j < len; j++) {
				k = ((k << 2) | base_bit_table[(int)seq[j]]) & kmask;
				if (j + 1 < ksize) continue;
				if (quality[j] <= '#' && quality[j+1] <= '#') continue; // trim ends with low qualities
				r = dna_rev_seq(k, ksize);
				if (r < k) {
					KMER.kmer = r;
				} else {
					KMER.kmer = k;
				}
				ret = get_kmerhash(hash, KMER);
				if (ret == NULL || ret->cnt < mincnt)  {
					continue;
				} 
				for (ii = 0; ii < ksize; ii++) {
					kmer_seq[ii] = bit_base_table[(ret->kmer >> ((ksize-1-ii) << 1)) & 0x03];
				}
				kmer_seq[ii] = '\0';
				if (ret->cnt2 < maxcnt2 && (float)ret->cnt2/(ret->cnt+ret->cnt2) < 0.01) { //TODO magic number
					pair = next_ref_pairv(pairs);
					pair->mt = SOMATIC;
					pair->pid = rid;
					/*
					pair->r1.name = strdup(read[0]->name.string);
					len_head = strlen(read[0]->header.string); 
					pair->r1.header = malloc(ksize + len_head + 2);
					memcpy(pair->r1.header, kmer_seq, ksize);
					pair->r1.header[ksize] = '_';
					memcpy(pair->r1.header+ksize+1, read[0]->header.string, len_head+1);
					//pair->r1.header = strdup(read[0]->header.string);
					pair->r1.seq = strdup(read[0]->seq.string);
					pair->r1.qual = strdup(read[0]->qual.string);
					pair->r2.name = strdup(read[1]->name.string);
					len_head = strlen(read[1]->header.string); 
					pair->r2.header = malloc(ksize + len_head + 2);
					memcpy(pair->r2.header, kmer_seq, ksize);
					pair->r2.header[ksize] = '_';
					memcpy(pair->r2.header+ksize+1, read[1]->header.string, len_head+1);
					//pair->r2.header = strdup(read[1]->header.string);
					pair->r2.seq = strdup(read[1]->seq.string);
					pair->r2.qual = strdup(read[1]->qual.string);
					*/
					//fprintf(stderr, "type=%d\t@%s\n%s\n+\n%s\n", pair->mt, pair->r1.header, pair->r1.seq, pair->r1.qual);
					//fprintf(stderr, "@%s\n%s\n+\n%s\n", pair->r2.header, pair->r2.seq, pair->r2.qual);
				}
				else if (ret->cnt2 > maxcnt2 * 2) { //TODO magic number
					if (is_somatic) continue;
					pair = next_ref_pairv(pairs);
					pair->mt = GERMLINE;
					pair->pid = rid;
					/*
					pair->r1.name = strdup(read[0]->name.string);
					len_head = strlen(read[0]->header.string); 
					pair->r1.header = malloc(ksize + len_head + 2);
					memcpy(pair->r1.header, kmer_seq, ksize);
					pair->r1.header[ksize] = '_';
					memcpy(pair->r1.header+ksize+1, read[0]->header.string, len_head+1);
					//pair->r1.header = strdup(read[0]->header.string);
					pair->r1.seq = strdup(read[0]->seq.string);
					pair->r1.qual = strdup(read[0]->qual.string);
					pair->r2.name = strdup(read[1]->name.string);
					len_head = strlen(read[1]->header.string); 
					pair->r2.header = malloc(ksize + len_head + 2);
					memcpy(pair->r2.header, kmer_seq, ksize);
					pair->r2.header[ksize] = '_';
					memcpy(pair->r2.header+ksize+1, read[1]->header.string, len_head+1);
					//pair->r2.header = strdup(read[1]->header.string);
					pair->r2.seq = strdup(read[1]->seq.string);
					pair->r2.qual = strdup(read[1]->qual.string);
					*/
				}
				//goto label;
			}
		}
	}
	if (read[1]) free_sequence(read[1]);
	return pairs;
}

void destroy_pairv(pairv *pairs) {
	//pair_t *pair;
	//uint32_t i;
/*
	for (i = 0; i < count_pairv(pairs); i++) {
		pair = ref_pairv(pairs, i);
		destroy_pair(&pair->r1);
		destroy_pair(&pair->r2);
	}
*/
	free_pairv(pairs);
}

static inline int cmp_pair_func(pair_t p1, pair_t p2, void *obj) {

	if (p1.pid == p2.pid) {
		if (p1.mt == p2.mt)
			return 0;
		else 
			return p1.mt-p2.mt;
	}


	return p1.pid-p2.pid;
	obj = obj;
}

define_quick_sort(sort_pairs, pair_t, cmp_pair_func);

inline int trim_lowq(char *qual, int len, char min) {
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

void dedup_pairs(pairv *pairs, FILE *out1, FILE *out2, FILE *out3, FILE *out4, FileReader *f1, FileReader *f2) {
	pair_t *pair = NULL;
	uint32_t i, rid = 0;
	//char pre1[256] = "", pre2[256] = "";
	//mut_type pre_t = SOMATIC;
	u32hash *somaids;
	u32hash *germids;
	uint32_t *skey, *gkey;
	int exists;
	Sequence *read[2];
	int is_fq[2];

	read[0] = read[1] = NULL;
	is_fq[0] = (guess_seq_file_type(f1) == 2);
	is_fq[1] = (guess_seq_file_type(f2) == 2);

	somaids = init_u32hash(1047);
	germids = init_u32hash(1047);

//	sort_pairs(ref_pairv(pairs, 0), count_pairv(pairs), NULL);
	
	for (i = 0; i < count_pairv(pairs); i++) {
		pair = ref_pairv(pairs, i);
		if (pair->mt == SOMATIC) {
			skey = prepare_u32hash(somaids, pair->pid, &exists);
			if (!exists) *skey = pair->pid;
		} else if (pair->mt == GERMLINE) {
			gkey = prepare_u32hash(germids, pair->pid, &exists);
			if (!exists) *gkey = pair->pid;
		}
	}
	while (1) {
		//label:
		if (!(is_fq[0]?fread_fastq(&read[0], f1):fread_fasta(&read[0], f1)))  break;
		if (!(is_fq[1]?fread_fastq(&read[1], f2):fread_fasta(&read[1], f2)))  break;
		rid ++;
		if ((rid & 0xFFFFU) == 0) {
			fprintf(stdout, "[%s] parsed %10u pairs\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		if (exists_u32hash(somaids, rid)) {
			i = trim_lowq(read[0]->qual.string, read[0]->qual.size, '#');
			if (i == 0) {
				fprintf(out1, "@%s\n%s\n+\n%s\n", read[0]->header.string, "N", "#"); 
			} else {
				fprintf(out1, "@%s\n%.*s\n+\n%.*s\n", read[0]->header.string, i, read[0]->seq.string, i, read[0]->qual.string); 
			}
			i = trim_lowq(read[1]->qual.string, read[1]->qual.size, '#');
			if (i == 0) {
				fprintf(out2, "@%s\n%s\n+\n%s\n", read[1]->header.string, "N", "#"); 
			} else {
				fprintf(out2, "@%s\n%.*s\n+\n%.*s\n", read[1]->header.string, i, read[1]->seq.string, i, read[1]->qual.string); 
			}
	//		fprintf(out2, "@%s\n%s\n+\n%s\n", read[1]->header.string, read[1]->seq.string, read[1]->qual.string);
		}
		if (exists_u32hash(germids, rid)) {
			i = trim_lowq(read[0]->qual.string, read[0]->qual.size, '#');
			if (i == 0) {
				fprintf(out3, "@%s\n%s\n+\n%s\n", read[0]->header.string, "N", "#"); 
			} else {
				fprintf(out3, "@%s\n%.*s\n+\n%.*s\n", read[0]->header.string, i, read[0]->seq.string, i, read[0]->qual.string); 
			}
			i = trim_lowq(read[1]->qual.string, read[1]->qual.size, '#');
			if (i == 0) {
				fprintf(out4, "@%s\n%s\n+\n%s\n", read[1]->header.string, "N", "#"); 
			} else {
				fprintf(out4, "@%s\n%.*s\n+\n%.*s\n", read[1]->header.string, i, read[1]->seq.string, i, read[1]->qual.string); 
			}
			//fprintf(out3, "@%s\n%s\n+\n%s\n", read[0]->header.string, read[0]->seq.string, read[0]->qual.string);
			//fprintf(out4, "@%s\n%s\n+\n%s\n", read[1]->header.string, read[1]->seq.string, read[1]->qual.string);
		}
	}
		/*
		if (pre1[0] == '\0') {
			memcpy(pre1, pair->r1.seq, strlen(pair->r1.seq)+1);
			memcpy(pre2, pair->r2.seq, strlen(pair->r2.seq)+1);
			pre_t = pair->mt;
			if (pair->mt == SOMATIC) {
				fprintf(out1, "@%s\n%s\n+\n%s\n", pair->r1.header, pair->r1.seq, pair->r1.qual);
				fprintf(out2, "@%s\n%s\n+\n%s\n", pair->r2.header, pair->r2.seq, pair->r2.qual);
			} else if (pair->mt == GERMLINE) {
				fprintf(out3, "@%s\n%s\n+\n%s\n", pair->r1.header, pair->r1.seq, pair->r1.qual);
				fprintf(out4, "@%s\n%s\n+\n%s\n", pair->r2.header, pair->r2.seq, pair->r2.qual);
			}
		} else {
			if (strcmp(pre1, pair->r1.seq) == 0 && strcmp(pre2, pair->r2.seq) == 0 && pre_t == pair->mt) {
				dups ++;
				continue;
			}
			else {
				if (pair->mt == SOMATIC) {
					fprintf(out1, "@%s\n%s\n+\n%s\n", pair->r1.header, pair->r1.seq, pair->r1.qual);
					fprintf(out2, "@%s\n%s\n+\n%s\n", pair->r2.header, pair->r2.seq, pair->r2.qual);
				} else if (pair->mt == GERMLINE) {
					fprintf(out3, "@%s\n%s\n+\n%s\n", pair->r1.header, pair->r1.seq, pair->r1.qual);
					fprintf(out4, "@%s\n%s\n+\n%s\n", pair->r2.header, pair->r2.seq, pair->r2.qual);
				}
				memcpy(pre1, pair->r1.seq, strlen(pair->r1.seq)+1);
				memcpy(pre2, pair->r2.seq, strlen(pair->r2.seq)+1);
				pre_t = pair->mt;
			}
		}
		*/
		//fprintf(stderr, "type%d\t%s\t%s\n", pair->mt, pair->r1.seq, pair->r2.seq);
//	printf("%s\t%s\n", pre1, pre2);
//	fprintf(stdout, "[%s] removed %u duplicated read pairs\n", __FUNCTION__, dups);
	if (read[1]) free_sequence(read[1]);
	free_u32hash(somaids);
	free_u32hash(germids);
	//fflush(stdout);
}

void cal_ctrl_kmers(kmerhash *hash, FileReader *fr, uint32_t ksize, int is_fq) {
	kmer_t KMER, *ret;
	Sequence *seq;
	uint64_t k, r, kmask;
	uint32_t i, len, rid;
//	kmask = (1LLU << (2 * ksize)) - 1;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	ret = 0;
	rid = 0;
	seq = NULL;
	while (is_fq?fread_fastq_adv(&seq, fr, 5):fread_fastq_adv(&seq, fr, 1)) {
		rid ++;
		if ((rid & 0xFFFFU) == 0) {
			fprintf(stdout, "[%s] parsed %10u reads\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		k = 0;
		len = seq->seq.size;
		for (i = 0; i < ksize-1; i++)
			k = (k << 2) | base_bit_table[(int)seq->seq.string[i]];
		for (i = 0; i <= len-ksize; i++) {
			k = ((k << 2) | base_bit_table[(int)seq->seq.string[i+ksize-1]])  & kmask;
			//if (i + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			ret = get_kmerhash(hash, KMER);
			if (ret) {
				if (ret->cnt2 < UNIQ_KMER_MAX_CNT)
					ret->cnt2 ++;
			}
		}
	}
	fprintf(stdout, "[%s] processed %10u reads\n", __FUNCTION__, rid);
	fflush(stdout);

	//return ret;
}
define_list(flist, char*);

inline int cmp_kmer(const void *e1, const void *e2) {
	kmer_t *k1, *k2;
	k1 = (kmer_t*)e1;
	k2 = (kmer_t*)e2;

	if (k1->cnt > k2->cnt) return -1;
	if (k1->cnt < k2->cnt) return 1;
	return 0;
}


int usage() {
	printf("clinsek - a tool for diagnosing known variations and discovering new (somatic) ones\n"
		   "Auther: Zechen Chong <chongzechen@gmail.com> & Wanding Zhou <zhouwanding@gmail.com>\n"
		   "Version: 1.01 (r20131023)\n"
		   "Usage:\n"
		   "  clinsek -1 <tumor_1.fq(.gz)> -2 <tumor_2.fq(.gz)> -3 <normal_1.fq(.gz)> -4 <normal_2.fq(.gz)> -r <reference> -o <output.kmer> [options]\n"
		   "Options:\n"
		   "  -h             This help\n"
		   "  -1 <string>    Treatment pair1 file in fastq format. Multiple treatment files could be input as -1 s1_1.fq -1 s2_1.fq ...\n"
		   "  -2 <string>    Treatment pair2 file in fastq format. Multiple treatment files could be input as -2 s1_2.fq -2 s2_1.fq ...\n"
		   "  -3 <string>    Control pair1 file in fastq format. Multiple treatment files could be input as -3 c1_1.fq -3 c2_1.fq ...\n"
		   "  -4 <string>    Control pair2 file in fastq format. Multiple treatment files could be input as -4 c1_2.fq -4 c2_2.fq ...\n"
		   "  -r <string>    Reference file in fasta format\n"
           "  -k <int>       Kmer size, <=31 [27]\n"
		   "  -o <string>    Output kmer\n"
		   "  -g <int>       Output germline events [0]\n"  
		   "  -m <int>       Minimum kmer count regarded as novo kmers [4]\n"
		   );

	return 1;
}

int main(int argc, char **argv) {
	kmerhash *khash;
	kmer_t KMER;
	pairv *pairs;
	FileReader *inf1, *inf2, *ctrlf1, *ctrlf2, *reff;
	FILE *out, *out1, *out2, *out3, *out4;
	char *in1file, *in2file, *outfile, *ctrl1file, *ctrl2file, *reffile;
	flist *in1list, *in2list, *ctrl1list, *ctrl2list;
	in1list = init_flist(2);
	in2list = init_flist(2);
	ctrl1list = init_flist(2);
	ctrl2list = init_flist(2);
	int c, is_fq, is_somatic = 1;
	uint32_t ksize = 27, mincnt = 4, i, maxcnt2 = 3;
	uint64_t ret;
	in1file = in2file = outfile = ctrl1file = ctrl2file = reffile = NULL;
	
	while ((c = getopt(argc, argv, "h1:2:3:4:k:o:m:r:g")) != -1) {
		switch (c) {
			case 'h': return usage();
			case '1': push_flist(in1list, optarg); break;
			case '2': push_flist(in2list, optarg); break;
			case '3': push_flist(ctrl1list, optarg); break;
			case '4': push_flist(ctrl2list, optarg); break;
			case 'k': ksize = atoi(optarg); break;
			case 'o': outfile = optarg; break;
			case 'm': mincnt = atoi(optarg); break;
			case 'r': reffile = optarg; break;
			case 'g': is_somatic = 0; break;
			default: return usage();
		}
	}
	if (count_flist(in1list) == 0 || count_flist(ctrl1list) == 0 || reffile == NULL) return usage();
	if (count_flist(in2list) != 0 && count_flist(in1list) != count_flist(in2list)) return usage();
	is_fq = 0;
	if ((reff = fopen_filereader(reffile)) == NULL) {
		fprintf(stderr, " -- Cannot open reference file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	if ((inf1 = fopen_m_filereader(count_flist(in1list), as_array_flist(in1list))) == NULL) {
		fprintf(stderr, " -- Cannot open input file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	} else {
		is_fq = guess_seq_file_type(inf1);
		switch (is_fq) {
			case 1: is_fq = 0; break;
			case 2: is_fq = 1; break;
			default: fprintf(stderr, "unknown file type\n");
			abort();
		}
	}
	if ((inf2 = fopen_m_filereader(count_flist(in2list), as_array_flist(in2list))) == NULL) {
		fprintf(stderr, " -- Cannot open input file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	} else {
		is_fq = guess_seq_file_type(inf2);
		switch (is_fq) {
			case 1: is_fq = 0; break;
			case 2: is_fq = 1; break;
			default: fprintf(stderr, "unknown file type\n");
			abort();
		}
	}
	if ((ctrlf1 = fopen_m_filereader(count_flist(ctrl1list), as_array_flist(ctrl1list))) == NULL) {
		fprintf(stderr, " -- Cannot open input file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		abort();
	} else {
		is_fq = guess_seq_file_type(ctrlf1);
		switch (is_fq) {
			case 1: is_fq = 0; break;
			case 2: is_fq = 1; break;
			default: fprintf(stderr, "unknown file type\n");
			abort();
		}
	}
	if ((ctrlf2 = fopen_m_filereader(count_flist(ctrl2list), as_array_flist(ctrl2list))) == NULL) {
		fprintf(stderr, " -- Cannot open input file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		abort();
	} else {
		is_fq = guess_seq_file_type(ctrlf2);
		switch (is_fq) {
			case 1: is_fq = 0; break;
			case 2: is_fq = 1; break;
			default: fprintf(stderr, "unknown file type\n");
			abort();
		}
	}
	if ((out = fopen(outfile, "w")) == NULL) {
		fprintf(stderr, " -- Please provide output file -o [output] --\n");
		fflush(stderr);
		abort();
	}

	if ((out1 = fopen("somatic_novo_kmer_read1.fq", "w")) == NULL) {
		fprintf(stderr, " cannot open file to write\n");
		fflush(stderr);
		abort();
	}

	if ((out2 = fopen("somatic_novo_kmer_read2.fq", "w")) == NULL) {
		fprintf(stderr, " cannot open file to write\n");
		fflush(stderr);
		abort();
	}
	if ((out3 = fopen("germline_novo_kmer_read1.fq", "w")) == NULL) {
		fprintf(stderr, " cannot open file to write\n");
		fflush(stderr);
		abort();
	}
	if ((out4 = fopen("germline_novo_kmer_read2.fq", "w")) == NULL) {
		fprintf(stderr, " cannot open file to write\n");
		fflush(stderr);
		abort();
	}
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "Building kmer...\n");
	fflush(stdout);
	khash = init_kmerhash(1023);
	khash = build_kmerhash(inf1, ksize, is_fq, khash);
	khash = build_kmerhash(inf2, ksize, is_fq, khash);
	fprintf(stdout, " %llu kmers loaded\n\n", (unsigned long long)count_kmerhash(khash));
	fflush(stdout);
	
	fclose_filereader(inf1);
	fclose_filereader(inf2);

	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "Filtering kmers from reference...\n");
	fflush(stdout);
	ret = filter_ref_kmers(khash, reff, ksize);
	fprintf(stdout, "Filtered %llu kmer from reference\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "There are still %llu kmers left\n\n", (unsigned long long)count_kmerhash(khash));
	fflush(stdout);
	
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "Calculating kmers from control...\n");
	fflush(stdout);
	cal_ctrl_kmers(khash, ctrlf1, ksize, is_fq);
	//fprintf(stdout, "Filtered %llu kmer from control file 1\n", (unsigned long long)ret);
	cal_ctrl_kmers(khash, ctrlf2, ksize, is_fq);
	//fprintf(stdout, "Filtered %llu kmer from control file 2\n", (unsigned long long)ret);
	//fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	//fprintf(stdout, "There are still %llu kmers left\n\n", (unsigned long long)count_kmerhash(khash));
	fflush(stdout);

	ret = 0;
	inf1 = fopen_m_filereader(count_flist(in1list), as_array_flist(in1list));
	inf2 = fopen_m_filereader(count_flist(in2list), as_array_flist(in2list));
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "Loading sequences...\n");
	fflush(stdout);
	pairs =  loadkmerseq(khash, ksize, mincnt, maxcnt2, inf1, inf2, is_somatic); //TODO
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "There are %llu pairs loaded\n\n", (unsigned long long)count_pairv(pairs));
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "Remove duplicate sequences...\n");
	fflush(stdout);
	//reset_filereader(inf1);
	//reset_filereader(inf2);
	fclose_filereader(inf1);
	fclose_filereader(inf2);
	inf1 = fopen_m_filereader(count_flist(in1list), as_array_flist(in1list));
	inf2 = fopen_m_filereader(count_flist(in2list), as_array_flist(in2list));
	dedup_pairs(pairs, out1, out2, out3, out4, inf1, inf2);
	fprintf(stdout, "[%s]\n\n", date()); fflush(stdout);
	fflush(stdout);
	reset_iter_kmerhash(khash);
	while (iter_kmerhash(khash, &KMER)) {
		//if (KMER.cnt < mincnt || KMER.cnt2 > maxcnt2) continue;
		if (KMER.cnt < mincnt) continue;
		if (KMER.cnt2 < maxcnt2 && (float)KMER.cnt2/(KMER.cnt+KMER.cnt2) < 0.01) {
			ret ++;
			for (i = 0; i < ksize; i++) {
				fprintf(out, "%c", bit_base_table[(KMER.kmer >> ((ksize-1-i) << 1)) & 0x03]);
			}
			fprintf(out, "\t%llu\t%llu\t%s\n", (unsigned long long)KMER.cnt, (unsigned long long)KMER.cnt2, "SOMATIC");
		} else if (KMER.cnt2 > maxcnt2*2) {
			ret ++;
			for (i = 0; i < ksize; i++) {
				fprintf(out, "%c", bit_base_table[(KMER.kmer >> ((ksize-1-i) << 1)) & 0x03]);
			}
			fprintf(out, "\t%llu\t%llu\t%s\n", (unsigned long long)KMER.cnt, (unsigned long long)KMER.cnt2, "GERMLINE");
		}
	}
	
	fprintf(stdout, "%llu kmers passed the minimum frequency cutoff in tumor (%u) and maximum frequency cutoff in normal (%u)\n", (unsigned long long)ret, mincnt, maxcnt2);
	fflush(stdout);
	destroy_pairv(pairs);
	free_kmerhash(khash);
	free_flist(in1list);
	free_flist(in2list);
	free_flist(ctrl1list);
	free_flist(ctrl2list);
	fclose_filereader(inf1);
	fclose_filereader(inf2);
	fclose_filereader(ctrlf1);
	fclose_filereader(ctrlf2);
	fclose_filereader(reff);
	fclose(out);
	fclose(out1);
	fclose(out2);
	fclose(out3);
	fclose(out4);
	fprintf(stderr, "Program exit normally\n");
	return 0;
}
