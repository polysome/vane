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

#include "echelonform.h"
#include "newton_john.h"
#include "permutation.h"
#include "trsm.h"
#include "ple.h"

rci_t mzed_echelonize(mzed_t *A, int full) {
	if (A->finite_field->degree > A->nrows) {
		return mzed_echelonize_naive(A, full);
	} else if ((A->nrows * A->ncols * A->w <= 2*__M4RIE_PLE_CUTOFF)) {
		return mzed_echelonize_newton_john(A, full);
	} else {
		return mzed_echelonize_ple(A, full);
	}
}

rci_t mzd_slice_echelonize_ple(mzd_slice_t *A, int full) {
	mzp_t *P = mzp_init(A->nrows);
	mzp_t *Q = mzp_init(A->ncols);
	rci_t r;

	if(full) {
		r = mzd_slice_pluq(A, P, Q);

		mzd_slice_t *U = mzd_slice_init_window(A, 0, 0, r, r);
		const rci_t r_radix = m4ri_radix*(r/m4ri_radix);

		if(r_radix == r && r!=A->ncols) {
			mzd_slice_t *B = mzd_slice_init_window(A, 0, r, r, A->ncols);
			for(rci_t i = 0; i < r; ++i)
				mzd_slice_write_elem(U, i, i, 1);
			mzd_slice_trsm_upper_left(U, B);
			mzd_slice_free_window(B);
		} else if (r_radix != r && r!=A->ncols) {
			assert(r_radix < r);

			if(A->ncols > r_radix+m4ri_radix) {
				mzd_slice_t *B0  = mzd_slice_submatrix(NULL, A, 0, r_radix, r, r_radix+m4ri_radix);
				mzd_slice_t *B0w = mzd_slice_init_window(    A, 0, r_radix, r, r_radix+m4ri_radix);
				mzd_slice_t *B1  = mzd_slice_init_window(A, 0, r_radix+m4ri_radix, r, A->ncols);

				for(rci_t i = 0; i < r; ++i)
					mzd_slice_write_elem(U, i, i, 1);

				mzd_slice_trsm_upper_left(U, B0);
				mzd_slice_trsm_upper_left(U, B1);

				mzd_slice_copy(B0w, B0);
				mzd_slice_free(B0);
				mzd_slice_free_window(B0w);
				mzd_slice_free_window(B1);

			} else {
				mzd_slice_t *B = mzd_slice_submatrix(NULL, A, 0, r_radix, r, A->ncols);
				mzd_slice_t *Bw = mzd_slice_init_window(A, 0, r_radix, r, A->ncols);

				for(rci_t i = 0; i < r; ++i)
					mzd_slice_write_elem(U, i, i, 1);

				mzd_slice_trsm_upper_left(U, B);

				mzd_slice_copy(Bw, B);
				mzd_slice_free_window(Bw);
				mzd_slice_free(B);     
			}
		}

		mzd_slice_set_ui(U, 1);

		mzd_slice_free_window(U);

		if(r) {
			mzd_slice_t *A0 = mzd_slice_init_window(A, 0, 0, r, A->ncols);
			mzd_slice_apply_p_right(A0, Q);
			mzd_slice_free_window(A0);
		}

	} else {
		r = mzd_slice_ple(A, P, Q);

		for(rci_t i = 0; i < r; ++i) {
			for(int e=0; e < A->depth; e++) {
				for(rci_t j = 0; j <= i; j++) {
					int const length = MIN(m4ri_radix, i - j + 1);
					mzd_clear_bits(A->x[e], i, j, length);
				}
			}
			mzd_slice_write_elem(A, i, Q->values[i], 1);
		}
	}

	if(r != A->nrows) {
		mzd_slice_t *R = mzd_slice_init_window(A, r, 0, A->nrows, A->ncols);
		mzd_slice_set_ui(R, 0);
		mzd_slice_free_window(R);
	}

	mzp_free(P);
	mzp_free(Q);

	return r;
}


