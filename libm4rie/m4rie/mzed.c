/******************************************************************************
 *
 *            M4RIE: Linear Algebra over GF(2^e)
 *
 *    Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
 *
 *  Distributed under the terms of the GNU General Public License (GEL)
 *  version 2 or higher.
 *
 *    This code is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    General Public License for more details.
 *
 *  The full text of the GPL is available at:
 *
 *                  http://www.gnu.org/licenses/
 ******************************************************************************/

#include <stdlib.h>

#include "config.h"
#include "mzed.h"
#include "strassen.h"
#include "mzd_slice.h"
#include "conversion.h"

mzed_t *mzed_init(const gf2e* k, rci_t m, rci_t n) {
	mzed_t *A = (mzed_t *)m4ri_mm_malloc(sizeof(mzed_t));

	A->finite_field = k;
	A->w = gf2e_degree_to_w(A->finite_field);
	A->nrows = m;
	A->ncols = n;
	A->x = mzd_init(m, A->w*n);
	return A;
}

void mzed_free(mzed_t *A) {
	mzd_free(A->x);
	m4ri_mm_free(A);
}

void mzed_randomize(mzed_t *A) {
	unsigned int bitmask = (1<<A->finite_field->degree)-1;
	for(rci_t r=0; r<A->nrows; r++) {
		for(rci_t c=0; c<A->ncols; c++) {
			mzed_write_elem(A,r,c, arc4random()&bitmask);
		}
	}
}

mzed_t *mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	if (A->nrows != B->nrows || A->ncols != B->ncols || A->finite_field != B->finite_field) {
		m4ri_die("mzed_add: rows, columns and fields must match.\n");
	}
	if (C == NULL) {
		C = mzed_init(A->finite_field, A->nrows, A->ncols);
	} else if (C != A) {
		if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != A->ncols) {
			m4ri_die("mzed_add: rows and columns of returned matrix must match.\n");
		}
	}
	mzd_add(C->x, A->x, B->x);
	return C;
}

mzed_t *_mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	mzd_add(C->x, A->x, B->x);
	return C;
}


mzed_t *_mzed_mul_init(mzed_t *C, const mzed_t *A, const mzed_t *B, int clear) {
	if (A->ncols != B->nrows || A->finite_field != B->finite_field) {
		m4ri_die("mzed_mul: rows, columns and fields must match.\n");
	}
	if (C == NULL) {
		C = mzed_init(A->finite_field, A->nrows, B->ncols);
	} else {
		if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != B->ncols) {
			m4ri_die("mzed_mul: rows and columns of returned matrix must match.\n");
		}
		if (clear)
			mzed_set_ui(C,0);
	}
	return C;
}

mzed_t *mzed_mul(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	C = _mzed_mul_init(C,A,B, TRUE);
	_mzed_mul(C, A, B);
	return C;
}

mzed_t *mzed_addmul(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	C = _mzed_mul_init(C,A,B, FALSE);
	_mzed_addmul(C, A, B);
	return C;
}

mzed_t *_mzed_mul(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	if (A->nrows >= 512 && A->ncols >= 512 && B->ncols >= 512)
		return _mzed_addmul_karatsuba(C, A, B);

	const rci_t cutoff = _mzed_strassen_cutoff(C, A, B);
	return _mzed_mul_strassen(C, A, B, cutoff);
}

mzed_t *_mzed_addmul(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	if (A->nrows >= 512 && A->ncols >= 512 && B->ncols >= 512)
		return _mzed_addmul_karatsuba(C, A, B);

	const rci_t cutoff = _mzed_strassen_cutoff(C, A, B);
	return _mzed_addmul_strassen(C, A, B, cutoff);
}

mzed_t *mzed_mul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	C = _mzed_mul_init(C,A,B, TRUE);
	return _mzed_mul_naive(C, A, B);
}

mzed_t *mzed_addmul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	C = _mzed_mul_init(C,A,B, FALSE);
	return _mzed_mul_naive(C, A, B);
}

mzed_t *_mzed_mul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	const gf2e* ff = C->finite_field;
	for (rci_t i=0; i<C->nrows; ++i) {
		for (rci_t j=0; j<C->ncols; ++j) {
			for (rci_t k=0; k<A->ncols; ++k) {
				mzed_add_elem(C, i, j, ff->mul(ff, mzed_read_elem(A,i, k), mzed_read_elem(B, k, j)));
			}
		}
	}
	return C;
}


mzed_t *mzed_mul_scalar(mzed_t *C, const word a, const mzed_t *B) {
	/**
	 * The algorithm proceeds as follows:
	 */
	if(C == NULL)
		C = mzed_init(B->finite_field, B->nrows, B->ncols);

	const gf2e *ff = B->finite_field;

	/**
	 * 0) If a direct approach would need less lookups we use that.
	 */

	if(ff->degree > 8 || B->nrows*B->ncols < 1<<17) {
		mzed_copy(C, B);
		for(rci_t i=0; i<B->nrows; i++)
			mzed_rescale_row(C, i, 0, a);
		return C;
	}

	/**
	 * 1) We generate a lookup table of 16-bit wide entries
	 */
	const word mask_16 = (1<<16)-1;
	const word *mul = (const word*)gf2e_t16_init(B->finite_field, a);

	/**
	 * 2) We use that lookup table to do 4 lookups per word
	 */

	for(rci_t i=0; i<C->nrows; i++) {
		word *c_row = C->x->rows[i];
		const word *b_row = B->x->rows[i];
		for(wi_t j=0; j<C->x->width-1; j++) {
			const word tmp = b_row[j];
			const word a0 = tmp & mask_16;
			const word a1 = tmp>>16 & mask_16;
			const word a2 = tmp>>32 & mask_16;
			const word a3 = tmp>>48 & mask_16;
			c_row[j] = mul[a3]<<48 | mul[a2]<<32 | mul[a1]<<16 | mul[a0];      
		}
		/* deal with rest */
		const word tmp = b_row[B->x->width-1] & B->x->high_bitmask;
		const word a0 = tmp & mask_16;
		const word a1 = tmp>>16 & mask_16;
		const word a2 = tmp>>32 & mask_16;
		const word a3 = tmp>>48 & mask_16;
		c_row[C->x->width-1] &= ~B->x->high_bitmask;
		c_row[C->x->width-1] |= mul[a3]<<48 | mul[a2]<<32 | mul[a1]<<16 | mul[a0];
	}
	gf2e_t16_free((word*)mul);
	return C;
}


mzed_t *mzed_copy(mzed_t *A, const mzed_t *B) {
	if (A == B)
		return A;
	if (A == NULL)
		A = mzed_init(B->finite_field, B->nrows, B->ncols);
	if (A->finite_field != B->finite_field || A->nrows != B->nrows || A->ncols != B->ncols) {
		m4ri_die("mzed_copy: target matrix has wrong dimensions or base field.");
	}
	mzd_copy(A->x, B->x);
	return A;
}

rci_t mzed_echelonize_naive(mzed_t *A, int full) {
	rci_t start_row,r,c,i,elim_start;
	word x = 0;

	rci_t nr = A->nrows;
	rci_t nc = A->ncols;

	const gf2e *ff = A->finite_field;

	start_row = 0;

	for(c=0; c<nc; c++) {
		for(r=start_row; r<nr; r++) {
			x = mzed_read_elem(A, r, c);
			if (x) {
				mzed_rescale_row(A, r, c, gf2e_inv(ff, x));
				mzd_row_swap(A->x, r, start_row);
				if (full)
					elim_start = 0;
				else
					elim_start = start_row + 1;
				for(i=elim_start; i<nr; i++) {
					if (i==start_row) 
						continue;
					x = mzed_read_elem(A,i,c);
					if(!x) continue;
					/* clear row */
					mzed_add_multiple_of_row(A, i, A, start_row, x, c);
				}
				start_row++;
				break;
			}
		}
	}
	return start_row;
}

void mzed_set_ui(mzed_t *A, word value) {
	mzd_set_ui(A->x, 0);
	if(!value)
		return;
	for(rci_t i=0; i< MIN(A->ncols,A->nrows); i++) {
		mzed_write_elem(A, i, i, value);
	}
}

void mzed_print(const mzed_t *A) {
	char formatstr[10];
	int width = (A->w/4);
	if (A->w%4) 
		width += 1;
	sprintf(formatstr,"%%%dx",width);
	for (rci_t i=0; i < A->nrows; ++i) {
		printf("[");
		for (rci_t j=0; j < A->ncols; j++) {
			word tmp = mzed_read_elem(A,i,j);
			printf(formatstr,(int)tmp);
			if(j<A->ncols-1)
				printf(" ");
		}
		printf("]\n");
	}
}

void mzed_add_multiple_of_row(mzed_t *A, rci_t ar, const mzed_t *B, rci_t br, word x, rci_t start_col) {
	assert(A->ncols == B->ncols && A->finite_field == B->finite_field);
	assert(start_col < A->ncols);

	const gf2e *ff = A->finite_field;

	if (x == 0) {
		return;
	} else if(x == 1) {
		mzed_add_row(A, ar, B, br, start_col);
		return;
	}

	const rci_t start = A->w*start_col;
	const wi_t startblock = start/m4ri_radix;
	const word bitmask_end = A->x->high_bitmask;

	mzd_t *from_x = B->x;
	mzd_t *to_x = A->x;
	word *_f = from_x->rows[br];
	word *_t = to_x->rows[ar];
	wi_t j;

	register word __f = _f[startblock]>>(start%m4ri_radix);
	register word __t = _t[startblock];

	if(A->w == 2) {
		switch( (start/2) % 32) {
			case  0:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 0;  __f >>= 2;
			case  1:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 2;  __f >>= 2;
			case  2:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 4;  __f >>= 2;
			case  3:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 6;  __f >>= 2;
			case  4:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 8;  __f >>= 2;
			case  5:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<10;  __f >>= 2;
			case  6:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<12;  __f >>= 2;
			case  7:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<14;  __f >>= 2;
			case  8:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<16;  __f >>= 2;
			case  9:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<18;  __f >>= 2;
			case 10:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<20;  __f >>= 2;
			case 11:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<22;  __f >>= 2;
			case 12:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<24;  __f >>= 2;
			case 13:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<26;  __f >>= 2;
			case 14:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<28;  __f >>= 2;
			case 15:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<30;  __f >>= 2;
			case 16:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<32;  __f >>= 2;
			case 17:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<34;  __f >>= 2;
			case 18:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<36;  __f >>= 2;
			case 19:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<38;  __f >>= 2;
			case 20:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<40;  __f >>= 2;
			case 21:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<42;  __f >>= 2;
			case 22:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<44;  __f >>= 2;
			case 23:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<46;  __f >>= 2;
			case 24:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<48;  __f >>= 2;
			case 25:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<50;  __f >>= 2;
			case 26:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<52;  __f >>= 2;
			case 27:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<54;  __f >>= 2;
			case 28:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<56;  __f >>= 2;
			case 29:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<58;  __f >>= 2;
			case 30:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<60;  __f >>= 2;
			case 31:  __t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<62;  break;
			default: m4ri_die("impossible");
		}

		if(to_x->width-startblock == 1) {
			_t[startblock] &= ~bitmask_end;
			_t[startblock] ^= __t & bitmask_end;
			return;
		} else {
			_t[startblock] = __t;
		}

		for(j=startblock+1; j<to_x->width -1; j++) {
			__f = _f[j], __t = _t[j];
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 0;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 2;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 4;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 6;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<< 8;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<10;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<12;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<14;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<16;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<18;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<20;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<22;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<24;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<26;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<28;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<30;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<32;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<34;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<36;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<38;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<40;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<42;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<44;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<46;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<48;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<50;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<52;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<54;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<56;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<58;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<60;  __f >>= 2;
			__t ^= ff->mul(ff, x, __f & 0x0000000000000003ULL)<<62;
			_t[j] = __t;
		}

		switch(to_x->ncols % m4ri_radix) {
			case  0: _t[j] ^= ff->mul(ff, x, (_f[j] & 0xC000000000000000ULL)>>62)<<62;
			case 62: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x3000000000000000ULL)>>60)<<60;
			case 60: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0C00000000000000ULL)>>58)<<58;
			case 58: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0300000000000000ULL)>>56)<<56;
			case 56: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00C0000000000000ULL)>>54)<<54;
			case 54: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0030000000000000ULL)>>52)<<52;
			case 52: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000C000000000000ULL)>>50)<<50;
			case 50: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0003000000000000ULL)>>48)<<48;
			case 48: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000C00000000000ULL)>>46)<<46;
			case 46: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000300000000000ULL)>>44)<<44;
			case 44: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000C0000000000ULL)>>42)<<42;
			case 42: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000030000000000ULL)>>40)<<40;
			case 40: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000C000000000ULL)>>38)<<38;
			case 38: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000003000000000ULL)>>36)<<36;
			case 36: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000C00000000ULL)>>34)<<34;
			case 34: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000300000000ULL)>>32)<<32;
			case 32: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000C0000000ULL)>>30)<<30;
			case 30: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000030000000ULL)>>28)<<28;
			case 28: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000000C000000ULL)>>26)<<26;
			case 26: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000003000000ULL)>>24)<<24;
			case 24: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000C00000ULL)>>22)<<22;
			case 22: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000300000ULL)>>20)<<20;
			case 20: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000000C0000ULL)>>18)<<18;
			case 18: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000030000ULL)>>16)<<16;
			case 16: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000000000C000ULL)>>14)<<14;
			case 14: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000003000ULL)>>12)<<12;
			case 12: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000000C00ULL)>>10)<<10;
			case 10: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000000300ULL)>> 8)<< 8;
			case  8: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000000000C0ULL)>> 6)<< 6;
			case  6: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000000030ULL)>> 4)<< 4;
			case  4: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000000000000CULL)>> 2)<< 2;
			case  2: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000000003ULL)>> 0)<< 0;
		};

	} else if(A->w == 4) {
		switch( (start/4) % 16 ) {
			case  0: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<< 0;  __f >>= 4;
			case  1: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<< 4;  __f >>= 4;
			case  2: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<< 8;  __f >>= 4;
			case  3: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<12;  __f >>= 4;
			case  4: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<16;  __f >>= 4;
			case  5: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<20;  __f >>= 4;
			case  6: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<24;  __f >>= 4;
			case  7: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<28;  __f >>= 4;
			case  8: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<32;  __f >>= 4;
			case  9: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<36;  __f >>= 4;
			case 10: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<40;  __f >>= 4;
			case 11: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<44;  __f >>= 4;
			case 12: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<48;  __f >>= 4;
			case 13: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<52;  __f >>= 4;
			case 14: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<56;  __f >>= 4;
			case 15: __t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<60;  break;
			default: m4ri_die("impossible");
		}

		if(to_x->width-startblock == 1) {
			_t[startblock] &= ~bitmask_end;
			_t[startblock] ^= __t & bitmask_end;
			return;
		} else {
			_t[startblock] = __t;
		}

		for(j=startblock+1; j<to_x->width -1; j++) {
			__f = _f[j], __t = _t[j];
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<< 0;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<< 4;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<< 8;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<12;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<16;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<20;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<24;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<28;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<32;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<36;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<40;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<44;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<48;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<52;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<56;  __f >>= 4;
			__t ^= ff->mul(ff, x, __f & 0x000000000000000FULL)<<60;
			_t[j] = __t;
		}

		switch(to_x->ncols % m4ri_radix) {
			case  0: _t[j] ^= ff->mul(ff, x, (_f[j] & 0xF000000000000000ULL)>>60)<<60;
			case 60: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0F00000000000000ULL)>>56)<<56;
			case 56: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00F0000000000000ULL)>>52)<<52;
			case 52: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000F000000000000ULL)>>48)<<48;
			case 48: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000F00000000000ULL)>>44)<<44;
			case 44: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000F0000000000ULL)>>40)<<40;
			case 40: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000F000000000ULL)>>36)<<36;
			case 36: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000F00000000ULL)>>32)<<32;
			case 32: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000F0000000ULL)>>28)<<28;
			case 28: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000000F000000ULL)>>24)<<24;
			case 24: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000F00000ULL)>>20)<<20;
			case 20: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000000F0000ULL)>>16)<<16;
			case 16: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000000000F000ULL)>>12)<<12;
			case 12: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000000F00ULL)>> 8)<< 8;
			case  8: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000000000F0ULL)>> 4)<< 4;
			case  4: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000000000000FULL)>> 0)<< 0;
		};

	} else if (A->w == 8) {
		register word __t0 ,__t1, __f0, __f1;

		__f0 = _f[startblock]>>(start%m4ri_radix), __t0 = _t[startblock];
		switch( (start/8) % 8 ) {
			case 0: __t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<< 0; __f0 >>= 8;
			case 1: __t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<< 8; __f0 >>= 8;
			case 2: __t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<16; __f0 >>= 8;
			case 3: __t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<24; __f0 >>= 8;
			case 4: __t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<32; __f0 >>= 8;
			case 5: __t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<40; __f0 >>= 8;
			case 6: __t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<48; __f0 >>= 8;
			case 7: __t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<56; break;
			default: m4ri_die("impossible");
		}

		if(to_x->width-startblock == 1) {
			_t[startblock] &= ~bitmask_end;
			_t[startblock] ^= __t0 & bitmask_end;
			return;
		} else {
			_t[startblock] = __t0;
		}

		for(j=startblock+1; j+2 < to_x->width; j+=2) {
			__f0 = _f[j], __t0 = _t[j];
			__f1 = _f[j+1], __t1 = _t[j+1];
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<< 0; __f0 >>= 8;
			__t1 ^= ff->mul(ff, x, __f1 & 0x00000000000000FFULL)<< 0; __f1 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<< 8; __f0 >>= 8;
			__t1 ^= ff->mul(ff, x, __f1 & 0x00000000000000FFULL)<< 8; __f1 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<16; __f0 >>= 8;
			__t1 ^= ff->mul(ff, x, __f1 & 0x00000000000000FFULL)<<16; __f1 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<24; __f0 >>= 8;
			__t1 ^= ff->mul(ff, x, __f1 & 0x00000000000000FFULL)<<24; __f1 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<32; __f0 >>= 8;
			__t1 ^= ff->mul(ff, x, __f1 & 0x00000000000000FFULL)<<32; __f1 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<40; __f0 >>= 8;
			__t1 ^= ff->mul(ff, x, __f1 & 0x00000000000000FFULL)<<40; __f1 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<48; __f0 >>= 8;
			__t1 ^= ff->mul(ff, x, __f1 & 0x00000000000000FFULL)<<48; __f1 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<56;
			__t1 ^= ff->mul(ff, x, __f1 & 0x00000000000000FFULL)<<56;
			_t[j+0] = __t0;
			_t[j+1] = __t1;
		}

		for(; j < to_x->width-1; j++) {
			__f0 = _f[j], __t0 = _t[j];
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<< 0; __f0 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<< 8; __f0 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<16; __f0 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<24; __f0 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<32; __f0 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<40; __f0 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<48; __f0 >>= 8;
			__t0 ^= ff->mul(ff, x, __f0 & 0x00000000000000FFULL)<<56;
			_t[j] = __t0;
		}

		switch(to_x->ncols % m4ri_radix) {
			case  0: _t[j] ^= ff->mul(ff, x, (_f[j] & 0xFF00000000000000ULL)>>56)<<56;
			case 56: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00FF000000000000ULL)>>48)<<48;
			case 48: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000FF0000000000ULL)>>40)<<40;
			case 40: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000FF00000000ULL)>>32)<<32;
			case 32: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000FF000000ULL)>>24)<<24;
			case 24: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000000000FF0000ULL)>>16)<<16;
			case 16: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000000000FF00ULL)>> 8)<< 8;
			case  8: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000000000FFULL)>> 0)<< 0;
		};

	} else if (A->w == 16) {
		mzd_t *from_x = B->x;
		mzd_t *to_x = A->x;
		word *_f = from_x->rows[br];
		word *_t = to_x->rows[ar];
		size_t j;
		register word __t, __f;

		__f = _f[startblock]>>(start%m4ri_radix), __t = _t[startblock];
		switch( (start/16)%4 ) {
			case 0: __t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<< 0; __f >>= 16;
			case 1: __t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<16; __f >>= 16;
			case 2: __t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<32; __f >>= 16;
			case 3: __t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<48; break;
			default: m4ri_die("impossible");
		}
		if(to_x->width-startblock == 1) {
			_t[startblock] &= ~bitmask_end;
			_t[startblock] ^= __t & bitmask_end;
			return;
		} else {
			_t[startblock] = __t;
		}

		for(j=startblock+1; j+4<to_x->width; j+=4) {
			__f = _f[j], __t = _t[j];
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<< 0; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<16; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<32; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<48;
			_t[j] = __t;

			__f = _f[j+1], __t = _t[j+1];
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<< 0; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<16; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<32; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<48;
			_t[j+1] = __t;


			__f = _f[j+2], __t = _t[j+2];
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<< 0; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<16; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<32; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<48;
			_t[j+2] = __t;

			__f = _f[j+3], __t = _t[j+3];
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<< 0; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<16; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<32; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<48;
			_t[j+3] = __t;
		}
		for( ; j<to_x->width-1; j++) {
			__f = _f[j], __t = _t[j];
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<< 0; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<16; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<32; __f >>= 16;
			__t ^= ff->mul(ff, x, __f & 0x000000000000FFFFULL)<<48;
			_t[j] = __t;
		}

		switch(to_x->ncols % m4ri_radix) {
			case  0: _t[j] ^= ff->mul(ff, x, (_f[j] & 0xFFFF000000000000ULL)>>48)<<48;
			case 48: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x0000FFFF00000000ULL)>>32)<<32;
			case 32: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x00000000FFFF0000ULL)>>16)<<16;
			case 16: _t[j] ^= ff->mul(ff, x, (_f[j] & 0x000000000000FFFFULL)>> 0)<< 0;
		};

	}  else {
		for(rci_t j=start_col; j<B->ncols; j++) {
			mzed_add_elem(A, ar, j, ff->mul(ff, x, mzed_read_elem(B, br, j)));
		}
	}
}
