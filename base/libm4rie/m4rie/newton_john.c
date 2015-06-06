/******************************************************************************
 *
 *            M4RIE: Linear Algebra over GF(2^e)
 *
 *    Copyright (C) 2010,2011 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include "config.h"

#include <m4ri/misc.h>
#include <m4ri/mzd.h>
#include <m4ri/brilliantrussian.h>
#include <m4ri/xor.h>

#include "newton_john.h"
#include "trsm.h"
#include "ple.h"
#include "conversion.h"

njt_mzed_t *njt_mzed_init(const gf2e *ff, const rci_t ncols) {
	njt_mzed_t *T =  m4ri_mm_malloc(sizeof(njt_mzed_t));
	T->L = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(ff->degree), sizeof(rci_t));
	T->T = mzed_init(ff, __M4RI_TWOPOW(ff->degree), ncols);
	T->M = mzed_init(ff, ff->degree, ncols);
	return T;
}

void njt_mzed_free(njt_mzed_t *T) {
	mzed_free(T->M);
	mzed_free(T->T);
	m4ri_mm_free(T->L);
	m4ri_mm_free(T);
}

/**
 * Compute C[rc,i] = C[rc,i] + T0[r0,i] + ... + T3[r3,i] for 0 <= i < ncols
 *
 * \param C Matrix
 * \apram rc Row index
 * \param T0 Matrix
 * \param r0 Row index
 * \param T1 Matrix
 * \param r1 Row index
 * \param T2 Matrix
 * \param r2 Row index
 * \param T3 Matrix
 * \param r3 Row index
 */

static inline void mzed_combine4(mzed_t *C, rci_t rc, 
		mzed_t *T0, rci_t r0, mzed_t *T1, rci_t r1, mzed_t *T2, rci_t r2, mzed_t *T3, rci_t r3) {
	const word *t[4] = {T0->x->rows[r0], T1->x->rows[r1], T2->x->rows[r2], T3->x->rows[r3]};
	_mzd_combine_4(C->x->rows[rc], t, C->x->width);
}

/**
 * Compute C[rc,i] = C[rc,i] + T0[r0,i] + ... + T7[r7,i] for 0 <= i < ncols
 *
 * \param C Matrix
 * \apram rc Row index
 * \param T0 Matrix
 * \param r0 Row index
 * \param T1 Matrix
 * \param r1 Row index
 * \param T2 Matrix
 * \param r2 Row index
 * \param T3 Matrix
 * \param r3 Row index
 * \param T4 Matrix
 * \param r4 Row index
 * \param T5 Matrix
 * \param r5 Row index
 * \param T6 Matrix
 * \param r6 Row index
 * \param T7 Matrix
 * \param r7 Row index
 */

static inline void mzed_combine8(mzed_t *C, rci_t rc, 
		mzed_t *T0, rci_t r0, mzed_t *T1, rci_t r1, mzed_t *T2, rci_t r2, mzed_t *T3, rci_t r3,
		mzed_t *T4, rci_t r4, mzed_t *T5, rci_t r5, mzed_t *T6, rci_t r6, mzed_t *T7, rci_t r7) {
	const word *t[8] = {T0->x->rows[r0], T1->x->rows[r1], T2->x->rows[r2], T3->x->rows[r3], 
		T4->x->rows[r4], T5->x->rows[r5], T6->x->rows[r6], T7->x->rows[r7]};
	_mzd_combine_8(C->x->rows[rc], t, C->x->width);
}


/**
 * \brief Perform Gaussian reduction to reduced row echelon form on a
 * submatrix.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row endrow (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

rci_t _mzed_gauss_submatrix_full(mzed_t *A, const rci_t r, const rci_t c, const rci_t end_row, int k) {
	rci_t i,j,l;
	rci_t start_row = r;
	int found;
	word tmp;

	const gf2e *ff = A->finite_field;

	for (j=c; j<c+k; j++) {
		found = 0;
		for (i=start_row; i< end_row; i++) {
			/* first we need to clear the first columns */
			for (l=0; l<j-c; l++) {
				tmp = mzed_read_elem(A, i, c+l);
				if (tmp) mzed_add_multiple_of_row(A, i, A, r+l, tmp, c+l);
			}
			/* pivot? */
			const word x = mzed_read_elem(A, i, j);
			if (x) {
				mzed_rescale_row(A, i, j, gf2e_inv(ff, x));
				mzd_row_swap(A->x, i, start_row);

				/* clear above */
				for (l=r; l<start_row; l++) {
					tmp = mzed_read_elem(A, l, j);
					if (tmp) mzed_add_multiple_of_row(A, l, A, start_row, tmp, j);
				}
				start_row++;
				found = 1;
				break;
			}
		}
		if (found==0) {
			return j - c;
		}
	}
	return j - c;
}


njt_mzed_t *mzed_make_table(njt_mzed_t *T, const mzed_t *A, const rci_t r, const rci_t c) {
	assert(m4ri_radix > A->finite_field->degree);
	if (T == NULL)
		T = njt_mzed_init(A->finite_field, A->ncols);

	mzd_set_ui(T->M->x,0);

#if 0
	for(rci_t i=0; i< T->T->nrows; i+=2) {
		T->L[i] = i;
		mzed_add_multiple_of_row(T->T, i, A, r, A->finite_field->mul[i], c);

		T->L[i+1] = i+1;
		mzed_copy_row(T->T, i+1, T->T, i);
		mzed_add_row(T->T, i+1, A, r, c);
	}
#else  
	const int degree = A->finite_field->degree;
	const wi_t homeblock = A->w*c / m4ri_radix;
	const wi_t wide = T->M->x->width - homeblock;
	const word bitmask_end = T->M->x->high_bitmask;
	wi_t j;

	for(int i=0; i<degree; i++) {
		mzed_add_multiple_of_row(T->M, i, A, r, 1ULL<<i, c);
	}

	for(rci_t i=1; i < T->T->nrows; ++i) {
		word *ti = T->T->x->rows[i] + homeblock;
		word *ti1 = T->T->x->rows[i-1] + homeblock;

		const rci_t rowneeded = m4ri_codebook[degree]->inc[i - 1];
		const int id = m4ri_codebook[degree]->ord[i];
		T->L[id] = i;

		word *m = T->M->x->rows[rowneeded] + homeblock;

		/* there might still be stuff left over from the previous table creation,
		   here we assume that this is at most 8 * m4ri_radix bits away. */
		switch (homeblock) {
			case 0:
				break;
			default:
			case 8: *(ti-7) = 0;
			case 7: *(ti-6) = 0;
			case 6: *(ti-5) = 0;
			case 5: *(ti-4) = 0;
			case 4: *(ti-3) = 0;
			case 3: *(ti-2) = 0;
			case 2: *(ti-2) = 0;
			case 1: *(ti-1) = 0;
		}

		for(j = 0; j + 8 <= wide - 1; j += 8) {
			*ti++ = *m++ ^ *ti1++;
			*ti++ = *m++ ^ *ti1++;
			*ti++ = *m++ ^ *ti1++;
			*ti++ = *m++ ^ *ti1++;
			*ti++ = *m++ ^ *ti1++;
			*ti++ = *m++ ^ *ti1++;
			*ti++ = *m++ ^ *ti1++;
			*ti++ = *m++ ^ *ti1++;
		}
		switch(wide - j) {
			case 8:  *ti++ = *m++ ^ *ti1++;
			case 7:  *ti++ = *m++ ^ *ti1++;
			case 6:  *ti++ = *m++ ^ *ti1++;
			case 5:  *ti++ = *m++ ^ *ti1++;
			case 4:  *ti++ = *m++ ^ *ti1++;
			case 3:  *ti++ = *m++ ^ *ti1++;
			case 2:  *ti++ = *m++ ^ *ti1++;
			case 1:  *ti++ = (*m++ ^ *ti1++) & bitmask_end;
		}
	}
#endif

	return T;
}

rci_t mzed_echelonize_newton_john(mzed_t *A, int full) {
	const gf2e* ff = A->finite_field;

	rci_t r,c;

	rci_t k = ff->degree;

	/* cf. mzd_echelonize_m4ri */
	rci_t kk = (rci_t)m4ri_opt_k(A->x->nrows, A->x->ncols, 0);
	if (kk>=7)
		kk = 7;
	if ( (6*(1<<kk)*A->ncols / 8.0) > __M4RI_CPU_L2_CACHE / 2.0 )
		kk -= 1;
	kk = (6*kk)/k;

	/* enforcing bounds */
	if (kk == 0)
		kk = 1;
	else if (kk > 6)
		kk = 6;

	rci_t kbar = 0;

	njt_mzed_t *T0 = njt_mzed_init(ff, A->ncols);
	njt_mzed_t *T1 = njt_mzed_init(ff, A->ncols);
	njt_mzed_t *T2 = njt_mzed_init(ff, A->ncols);
	njt_mzed_t *T3 = njt_mzed_init(ff, A->ncols);
	njt_mzed_t *T4 = njt_mzed_init(ff, A->ncols);
	njt_mzed_t *T5 = njt_mzed_init(ff, A->ncols);

	r = 0;
	c = 0;
	while(c < A->ncols) {
		if(c+kk > A->ncols) kk = A->ncols - c;

		/**
		 * \todo we don't really compute the upper triangular form yet,
		 *       we need to implement _mzed_gauss_submatrix() and a better
		 *       table creation for that.
		 */
		kbar = _mzed_gauss_submatrix_full(A, r, c, A->nrows, kk);

		if (kbar == 6)  {
			mzed_make_table(T0, A, r,   c);
			mzed_make_table(T1, A, r+1, c+1);
			mzed_make_table(T2, A, r+2, c+2);
			mzed_make_table(T3, A, r+3, c+3);
			mzed_make_table(T4, A, r+4, c+4);
			mzed_make_table(T5, A, r+5, c+5);
			if(kbar == kk)
				mzed_process_rows6( A, r+6, A->nrows, c, T0, T1, T2, T3, T4, T5);
			if(full)
				mzed_process_rows6( A,   0,        r, c, T0, T1, T2, T3, T4, T5);
		} else if(kbar == 5) {
			mzed_make_table(T0, A, r,     c);
			mzed_make_table(T1, A, r+1, c+1);
			mzed_make_table(T2, A, r+2, c+2);
			mzed_make_table(T3, A, r+3, c+3);
			mzed_make_table(T4, A, r+4, c+4);
			if(kbar == kk)
				mzed_process_rows5( A, r+5, A->nrows, c, T0, T1, T2, T3, T4);
			if(full)
				mzed_process_rows5( A,   0,        r, c, T0, T1, T2, T3, T4);

		} else if(kbar == 4) {
			mzed_make_table(T0, A, r,   c);
			mzed_make_table(T1, A, r+1, c+1);
			mzed_make_table(T2, A, r+2, c+2);
			mzed_make_table(T3, A, r+3, c+3);
			if(kbar == kk)
				mzed_process_rows4( A, r+4, A->nrows, c, T0, T1, T2, T3);
			if(full)
				mzed_process_rows4( A,   0,        r, c, T0, T1, T2, T3);

		} else if(kbar == 3) {
			mzed_make_table(T0,  A, r,   c );
			mzed_make_table(T1,  A, r+1, c+1);
			mzed_make_table(T2,  A, r+2, c+2);
			if(kbar == kk)
				mzed_process_rows3( A, r+3, A->nrows, c, T0, T1, T2);
			if(full)
				mzed_process_rows3( A,   0,        r, c, T0, T1, T2);

		} else if(kbar == 2) {
			mzed_make_table(T0, A, r,   c );
			mzed_make_table(T1, A, r+1, c+1);
			if(kbar == kk)
				mzed_process_rows2( A, r+2, A->nrows, c, T0, T1);
			if(full)
				mzed_process_rows2( A,   0,        r, c, T0, T1);

		} else if (kbar == 1) {
			mzed_make_table(T0, A, r, c);
			if(kbar == kk)
				mzed_process_rows( A, r+1, A->nrows, c, T0);
			if(full)
				mzed_process_rows( A,   0,        r, c, T0);

		} else {
			c++;
		}
		r += kbar;
		c += kbar;
	}

	njt_mzed_free(T0);
	njt_mzed_free(T1);
	njt_mzed_free(T2);
	njt_mzed_free(T3);
	njt_mzed_free(T4);
	njt_mzed_free(T5);
	return r;
}

rci_t mzed_ple_newton_john(mzed_t *A, mzp_t *P, mzp_t *Q) {
	rci_t col_pos = 0;
	rci_t row_pos = 0;
	word tmp = 0;
	const gf2e *ff = A->finite_field;
	rci_t i,j;
	int found = 0;

	njt_mzed_t *T0 = njt_mzed_init(A->finite_field, A->ncols);

	while (row_pos < A->nrows && col_pos < A->ncols) {
		found = 0;
		for(j=col_pos; j<A->ncols; j++) {
			for(i=row_pos; i<A->nrows; i++) {
				if( (tmp = mzed_read_elem(A, i,j)) != 0) {
					found = 1;
					break;
				}
			}
			if (found)
				break;
		}
		if (found) {
			P->values[row_pos] = i;
			Q->values[row_pos] = j;
			mzed_row_swap(A, row_pos, i);

			if (j+1 < A->ncols) {
				mzed_rescale_row(A, row_pos, j+1, gf2e_inv(ff, tmp));
				mzed_make_table(T0, A, row_pos, j+1);
				mzed_process_rows(A, row_pos+1, A->nrows, j, T0);
			}
			row_pos++;
			col_pos = j + 1;
		} else {
			break;
		}
	}
	for (rci_t i = row_pos; i < A->nrows; ++i)
		P->values[i] = i;
	for (rci_t i = row_pos; i < A->ncols; ++i)
		Q->values[i] = i;
	for (rci_t i=0; i < row_pos; i++) {
		mzed_col_swap_in_rows(A, i, Q->values[i], i, A->nrows);
	}
	njt_mzed_free(T0);

	return row_pos;
}

mzed_t *_mzed_mul_newton_john0(mzed_t *C, const mzed_t *A, const mzed_t *B) {

	njt_mzed_t *T0 = njt_mzed_init(B->finite_field, B->ncols);

	for(rci_t i=0; i < A->ncols; i++) {
		mzed_make_table(T0, B, i, 0);
		for(rci_t j=0; j<A->nrows; j++)
			mzd_combine(C->x, j, 0, C->x, j, 0, T0->T->x, T0->L[mzed_read_elem(A, j, i)], 0);
	}
	njt_mzed_free(T0);
	return C;
}

mzed_t *_mzed_mul_newton_john(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	if (A->finite_field->degree > A->nrows)
		return _mzed_mul_naive(C, A, B);

	njt_mzed_t *T0 = njt_mzed_init(B->finite_field, B->ncols);
	njt_mzed_t *T1 = njt_mzed_init(B->finite_field, B->ncols);
	njt_mzed_t *T2 = njt_mzed_init(B->finite_field, B->ncols);
	njt_mzed_t *T3 = njt_mzed_init(B->finite_field, B->ncols);
	njt_mzed_t *T4 = njt_mzed_init(B->finite_field, B->ncols);
	njt_mzed_t *T5 = njt_mzed_init(B->finite_field, B->ncols);
	njt_mzed_t *T6 = njt_mzed_init(B->finite_field, B->ncols);
	njt_mzed_t *T7 = njt_mzed_init(B->finite_field, B->ncols);

	const rci_t kk = 8;
	const rci_t end = A->ncols/kk;

	rci_t blocksize = 1ULL<<30;

	if (A->nrows >= A->w*__M4RI_MUL_BLOCKSIZE) 
		blocksize = __M4RI_MUL_BLOCKSIZE/A->w; 

	rci_t giantstep, babystep;

	for (giantstep=0; giantstep + blocksize <= A->nrows; giantstep += blocksize) {
		for(rci_t i=0; i < end; i++) {
			mzed_make_table(T0, B, kk*i  , 0);
			mzed_make_table(T1, B, kk*i+1, 0);
			mzed_make_table(T2, B, kk*i+2, 0);
			mzed_make_table(T3, B, kk*i+3, 0);
			mzed_make_table(T4, B, kk*i+4, 0);
			mzed_make_table(T5, B, kk*i+5, 0);
			mzed_make_table(T6, B, kk*i+6, 0);
			mzed_make_table(T7, B, kk*i+7, 0);
			for(babystep = 0; babystep < blocksize; babystep++) {
				const rci_t j = giantstep + babystep;
				const rci_t x0 = T0->L[mzed_read_elem(A, j, kk*  i)];
				const rci_t x1 = T1->L[mzed_read_elem(A, j, kk*i+1)];
				const rci_t x2 = T2->L[mzed_read_elem(A, j, kk*i+2)];
				const rci_t x3 = T3->L[mzed_read_elem(A, j, kk*i+3)];
				const rci_t x4 = T4->L[mzed_read_elem(A, j, kk*i+4)];
				const rci_t x5 = T5->L[mzed_read_elem(A, j, kk*i+5)];
				const rci_t x6 = T6->L[mzed_read_elem(A, j, kk*i+6)];
				const rci_t x7 = T7->L[mzed_read_elem(A, j, kk*i+7)];
				mzed_combine8(C, j, T0->T, x0, T1->T, x1, T2->T, x2, T3->T, x3, T4->T, x4, T5->T, x5, T6->T, x6, T7->T, x7);
			}
		}
	}

	/* last giant step */
	for(rci_t i=0; i < end; i++) {
		mzed_make_table(T0, B, kk*i  , 0);
		mzed_make_table(T1, B, kk*i+1, 0);
		mzed_make_table(T2, B, kk*i+2, 0);
		mzed_make_table(T3, B, kk*i+3, 0);
		mzed_make_table(T4, B, kk*i+4, 0);
		mzed_make_table(T5, B, kk*i+5, 0);
		mzed_make_table(T6, B, kk*i+6, 0);
		mzed_make_table(T7, B, kk*i+7, 0);
		for(babystep = 0; babystep < A->nrows - giantstep; babystep++) {
			const rci_t j = giantstep + babystep;
			const rci_t x0 = T0->L[mzed_read_elem(A, j, kk*  i)];
			const rci_t x1 = T1->L[mzed_read_elem(A, j, kk*i+1)];
			const rci_t x2 = T2->L[mzed_read_elem(A, j, kk*i+2)];
			const rci_t x3 = T3->L[mzed_read_elem(A, j, kk*i+3)];
			const rci_t x4 = T4->L[mzed_read_elem(A, j, kk*i+4)];
			const rci_t x5 = T5->L[mzed_read_elem(A, j, kk*i+5)];
			const rci_t x6 = T6->L[mzed_read_elem(A, j, kk*i+6)];
			const rci_t x7 = T7->L[mzed_read_elem(A, j, kk*i+7)];
			mzed_combine8(C, j, T0->T, x0, T1->T, x1, T2->T, x2, T3->T, x3, T4->T, x4, T5->T, x5, T6->T, x6, T7->T, x7);
		}
	}

	if (A->ncols%kk) {
		for(rci_t i=kk*end; i < A->ncols; i++) {
			mzed_make_table(T0, B, i, 0);
			for(rci_t j=0; j<A->nrows; j++)
				mzd_combine(C->x, j, 0, C->x, j, 0, T0->T->x, T0->L[mzed_read_elem(A, j, i)], 0);
		}
	}

	njt_mzed_free(T0); njt_mzed_free(T1);  njt_mzed_free(T2); njt_mzed_free(T3);
	njt_mzed_free(T4); njt_mzed_free(T5);  njt_mzed_free(T6); njt_mzed_free(T7);
	return C;
}


mzed_t *mzed_mul_newton_john(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	C = _mzed_mul_init(C,A,B, TRUE);
	return _mzed_mul_newton_john(C, A, B);
}

mzed_t *mzed_addmul_newton_john(mzed_t *C, const mzed_t *A, const mzed_t *B) {
	C = _mzed_mul_init(C,A,B, FALSE);
	return _mzed_mul_newton_john(C, A, B);
}

mzed_t *mzed_invert_newton_john(mzed_t *B, const mzed_t *A) {
	assert(A->nrows == A->ncols);
	mzed_t *I = mzed_init(A->finite_field, A->nrows, A->ncols);
	mzed_set_ui(I, 1);
	mzed_t *T = mzed_concat(NULL, A, I);
	mzed_free(I);

	rci_t r = mzed_echelonize_newton_john(T, 1);
	if (r != A->nrows) 
		m4ri_die("mzed_invert_newton_john: input matrix does not have full rank.");
	B = mzed_submatrix(B, T, 0, A->ncols, A->nrows, T->ncols);
	mzed_free(T);
	return B;
}

void mzed_trsm_lower_left_newton_john(const mzed_t *L, mzed_t *B) {
	assert(L->finite_field == B->finite_field);
	assert(L->nrows == L->ncols);
	assert(B->nrows == L->ncols);

	const gf2e *ff = L->finite_field;
	if (__M4RI_TWOPOW(ff->degree) >= L->nrows) {
		mzed_trsm_lower_left_naive(L, B);
		return;
	}

	njt_mzed_t *T0 = njt_mzed_init(B->finite_field, B->ncols);

	for(rci_t i=0; i<B->nrows; i++) {
		mzed_rescale_row(B, i, 0, gf2e_inv(ff, mzed_read_elem(L, i, i)));
		mzed_make_table(T0, B, i, 0);
		for(rci_t j=i+1; j<B->nrows; j++)
			mzd_combine(B->x, j, 0, B->x, j, 0, T0->T->x, T0->L[mzed_read_elem(L, j, i)], 0);
	}
	njt_mzed_free(T0);
}

void mzed_trsm_upper_left_newton_john(const mzed_t *U, mzed_t *B) {
	assert(U->finite_field == B->finite_field);
	assert(U->nrows == U->ncols);
	assert(B->nrows == U->ncols);

	const gf2e *ff = U->finite_field;
	if ( (__M4RI_TWOPOW(ff->degree) >= U->nrows) ) {
		mzed_trsm_upper_left_naive(U, B);
		return;
	}

	njt_mzed_t *T0 = njt_mzed_init(B->finite_field, B->ncols);

	for(int i=B->nrows-1; i>=0; i--) {
		mzed_rescale_row(B, i, 0, gf2e_inv(ff, mzed_read_elem(U, i, i)));
		mzed_make_table(T0, B, i, 0);
		for(rci_t j=0; j<i; j++)
			mzd_combine(B->x, j, 0, B->x, j, 0, T0->T->x, T0->L[mzed_read_elem(U, j, i)], 0);
	}
	njt_mzed_free(T0);
}

void mzd_slice_trsm_lower_left_newton_john(const mzd_slice_t *L, mzd_slice_t *B) {
	assert(L->finite_field == B->finite_field);
	assert(L->nrows == L->ncols);
	assert(B->nrows == L->ncols);

	const gf2e *ff = L->finite_field;
	if (__M4RI_TWOPOW(ff->degree) >= L->nrows) {
		mzd_slice_trsm_lower_left_naive(L, B);
		return;
	}

	mzed_t *Be = mzed_cling(NULL, B);
	njt_mzed_t *T0 = njt_mzed_init(B->finite_field, B->ncols);

	for(rci_t i=0; i<B->nrows; i++) {
		mzed_rescale_row(Be, i, 0, gf2e_inv(ff, mzd_slice_read_elem(L, i, i)));
		mzed_make_table(T0, Be, i, 0);
		for(rci_t j=i+1; j<Be->nrows; j++)
			mzd_combine(Be->x, j, 0, Be->x, j, 0, T0->T->x, T0->L[mzd_slice_read_elem(L, j, i)], 0);
	}
	mzed_slice(B, Be);
	mzed_free(Be);
	njt_mzed_free(T0);
}

void mzd_slice_trsm_upper_left_newton_john(const mzd_slice_t *U, mzd_slice_t *B) {
	assert(U->finite_field == B->finite_field);
	assert(U->nrows == U->ncols);
	assert(B->nrows == U->ncols);

	const gf2e *ff = U->finite_field;
	if ( (__M4RI_TWOPOW(ff->degree) >= U->nrows)) {
		mzd_slice_trsm_upper_left_naive(U, B);
		return;
	}

	mzed_t *Be = mzed_cling(NULL, B);
	njt_mzed_t *T0 = njt_mzed_init(Be->finite_field, Be->ncols);

	for(int i=B->nrows-1; i>=0; i--) {
		mzed_rescale_row(Be, i, 0, gf2e_inv(ff, mzd_slice_read_elem(U, i, i)));
		mzed_make_table(T0, Be, i, 0);
		for(rci_t j=0; j<i; j++)
			mzd_combine(Be->x, j, 0, Be->x, j, 0, T0->T->x, T0->L[mzd_slice_read_elem(U, j, i)], 0);
	}
	mzed_slice(B, Be);
	mzed_free(Be);
	njt_mzed_free(T0);
}
