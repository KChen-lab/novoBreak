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
 
#ifndef __BITS_VEC_RJ_H
#define __BITS_VEC_RJ_H

#include <stdint.h>
#include <string.h>
#include <stdlib.h>

/* Useful functions when n_bit > 8 */

static inline void u8byte2bits(uint64_t val, uint8_t *dat, uint64_t offset, uint8_t size){
	uint8_t i, *v;
	v = &val;
#if __BYTE_ORDER == 1234
	for(i=0;i<size;i++){
		if(v[i >> 3] & (1U << (i & 0x7U))) dat[offset >> 3] |= 1U << (offset & 0x7U);
		else dat[offset >> 3] &= ~(1U << (offset & 0x7U));
		offset ++;
	}
#else
	for(i=0;i<size;i++){
		if(v[7 - (i >> 3)] & (1U << (i & 0x7U))) dat[offset >> 3] |= 1U << (offset & 0x7U);
		else dat[offset >> 3] &= ~(1U << (offset & 0x7U));
		offset ++;
	}
#endif
}

static inline uint64_t bits2u8byte(uint8_t *dat, uint64_t offset, uint8_t size){
	uint64_t ret;
	uint8_t i, *v;
	ret = 0;
	v = (uint8_t*)&ret;
#if __BYTE_ORDER == 1234
	for(i=0;i<size;i++){
		if(dat[offset >> 3] & (1U << (offset & 0x7U))) v[i >> 3] |= 1U << (i & 0x7U);
		offset ++;
	}
#else
	for(i=0;i<size;i++){
		if(dat[offset >> 3] & (1U << (offset & 0x7U))) v[7 - (i >> 3)] |= 1U << (i & 0x7U);
		offset ++;
	}
#endif
	return ret;
}

typedef struct {
	uint8_t *bits;
	uint64_t size;
	uint64_t cap;
	uint32_t n_bit;
} BitsVec;

static inline BitsVec* init_bitsvec(uint64_t size, uint32_t n_bit){
	BitsVec *vec;
	if(n_bit == 0) n_bit = 1;
	if(size < 8) size = 8;
	vec = malloc(sizeof(BitsVec));
	vec->n_bit = n_bit;
	vec->size  = 0;
	vec->cap   = size;
	vec->bits  = calloc(1, vec->size * vec->n_bit / 8);
	return vec;
}

static inline void free_bitsvec(BitsVec *vec){
	free(vec->bits);
	free(vec);
}

static inline int encap_bitsvec(BitsVec *vec, uint32_t n){
	uint64_t cap;
	if(vec->size + n <= vec->cap) return 0;
	cap = vec->cap;
	while(vec->size + n > vec->cap){
		if(vec->cap < 1024 * 1024){
			vec->cap <<= 1;
		} else {
			vec->cap += 1024 * 1024;
		}
	}
	vec->bits = realloc(vec->bits, (vec->cap * vec->n_bit) / 8);
	memset(vec->bits + (cap * vec->n_bit) / 8, 0, ((vec->cap - cap) * vec->n_cap) / 8);
	return 1;
}

static inline void set_m_bitsvec(BitsVec *vec, uint64_t idx, uint64_t m, uint8_t *dat, uint64_t offset){
	uint64_t i, off;
	off = (idx * vec->n_bit);
	for(i=0;i<vec->n_bit * m;i++){
		if(dat[offset >> 3] & (1U << (offset & 0x7U))) vec->bits[off >> 3] |= 1U << (off & 0x7U);
		offset ++;
		off ++;
	}
}

#define set_bitsvec(vec, dat, offset) set_m_bitsvec(vec, 1, dat, offset)

static inline void push_m_bitsvec(BitsVec *vec, uint64_t m, uint8_t *dat, uint64_t offset){
	encap_bitsvec(vec, m);
	set_m_bitsvec(vec, vec->size, m, dat, offset);
	vec->size ++;
}

#define push_bitsvec(vec, dat, offset) push_m_bitsvec(vec, 1, dat, offset)

static inline void get_m_bitsvec(BitsVec *vec, uint64_t idx, uint64_t m, uint8_t dat, uint64_t offset){
	uint64_t i, off;
	off = (idx * vec->n_bit);
	for(i=0;i<vec->n_bit * m;i++){
		if(vec->bits[off >> 3] & (1U << (off & 0x7U))){
			dat[offset >> 3] |= 1U << (offset & 0x7U);
		} else {
			dat[offset >> 3] &= ~(1U << (offset & 0x7U));
		}
		offset ++;
		off ++;
	}
}

#define get_bitsvec(vec, dat, offset) get_m_bitsvec(vec, 1, dat, offset)

static inline int pop_m_bitvec(BitsVec *vec, uint64_t m, uint8_t *dat, uint64_t offset){
	if(vec->size < m) return 0;
	vec->size -= m;
	get_m_bitsvec(vec, vec->size, m, dat, offset);
	return 1;
}

#define pop_bitsvec(vec, dat, offset) pop_m_bitsvec(vec, 1, dat, offset)

static inline void append_bitsvec(BitsVec *dst, BitsVec *src){
	encap_bitsvec(dst, src->size);
	push_m_bitsvec(dst, src->size, src->bits, offset);
}

#endif
