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
	uint32_t i, len, rid;
	int exists;
	Sequence *seq;
	rid = 0;
	KMER.cnt = 0;
	KMER.kmer = 0;
//	kmask = (1LLU << (2 * ksize)) - 1;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	seq = NULL;
	while (is_fq?fread_fastq_adv(&seq, fr, 5):fread_fastq_adv(&seq, fr, 1)) {
		rid ++;
		if ((rid & 0xFFFFU) == 0) {
			fprintf(stdout, "[%s] parsed %10u reads\r", __FUNCTION__, rid);
			fflush(stdout);
		}
		k = 0;
		len = seq->seq.size;
		for (i = 0; i < ksize-1; i++) {
			k = (k << 2) | base_bit_table[(int)seq->seq.string[i]];
		}
		for (i = 0; i <= len-ksize; i++) {
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

uint64_t filter_ctrl_kmers(kmerhash *hash, FileReader *fr, uint32_t ksize, int is_fq) {
	kmer_t KMER;
	Sequence *seq;
	uint64_t k, r, kmask, ret;
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
			ret += remove_kmerhash(hash, KMER);
		}
	}
	fprintf(stdout, "[%s] processed %10u reads\n", __FUNCTION__, rid);
	fflush(stdout);

	return ret;
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
		   "  -m <int>       Minimum kmer count regarded as novo kmers [4]\n"
		   );

	return 1;
}

int main(int argc, char **argv) {
	kmerhash *khash;
	kmer_t KMER;
	FileReader *inf1, *inf2, *ctrlf1, *ctrlf2, *reff;
	FILE *out;
	char *in1file, *in2file, *outfile, *ctrl1file, *ctrl2file, *reffile;
	flist *in1list, *in2list, *ctrl1list, *ctrl2list;
	in1list = init_flist(2);
	in2list = init_flist(2);
	ctrl1list = init_flist(2);
	ctrl2list = init_flist(2);
	int c, is_fq;
	uint32_t ksize = 27, mincnt = 4, i;
	uint64_t ret;
	in1file = in2file = outfile = ctrl1file = ctrl2file = reffile = NULL;
	
	while ((c = getopt(argc, argv, "h1:2:3:4:k:o:m:r:")) != -1) {
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

	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "Building kmer...\n");
	fflush(stdout);
	khash = init_kmerhash(1023);
	khash = build_kmerhash(inf1, ksize, is_fq, khash);
	khash = build_kmerhash(inf2, ksize, is_fq, khash);
	fprintf(stdout, " %llu kmers loaded\n\n", (unsigned long long)count_kmerhash(khash));
	fflush(stdout);
	
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "Filtering kmers from control...\n");
	fflush(stdout);
	ret = filter_ctrl_kmers(khash, ctrlf1, ksize, is_fq);
	fprintf(stdout, "Filtered %llu kmer from control file 1\n", (unsigned long long)ret);
	ret = filter_ctrl_kmers(khash, ctrlf2, ksize, is_fq);
	fprintf(stdout, "Filtered %llu kmer from control file 2\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "There are still %llu kmers left\n\n", (unsigned long long)count_kmerhash(khash));
	fflush(stdout);

	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "Filtering kmers from reference...\n");
	fflush(stdout);
	ret = filter_ref_kmers(khash, reff, ksize);
	fprintf(stdout, "Filtered %llu kmer from reference\n", (unsigned long long)ret);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	fprintf(stdout, "There are still %llu kmers left\n\n", (unsigned long long)count_kmerhash(khash));
	fflush(stdout);
//	qsort(khash->array, (unsigned long long)count_kmerhash(khash), sizeof(kmer_t), cmp_kmer);
	ret = 0;
	reset_iter_kmerhash(khash);
	while (iter_kmerhash(khash, &KMER)) {
		if (KMER.cnt < mincnt) continue;
		ret ++;
		for (i = 0; i < ksize; i++) {
			fprintf(out, "%c", bit_base_table[(KMER.kmer >> ((ksize-1-i) << 1)) & 0x03]);
		}
//		bits2seq(seq, &KMER.kmer, 0, ksize);
		fprintf(out, "\t%llu\n", (unsigned long long)KMER.cnt);
	}
	
	fprintf(stdout, "%llu kmers passed the minimum frequency cutoff (%u)\n", (unsigned long long)ret, mincnt);
	fflush(stdout);
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
	fprintf(stderr, "Program exit normally\n");
	return 0;
}
