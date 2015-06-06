/******************************************************************************
 *
 *            M4RIE: Linear Algebra over GF(2^e)
 *
 *    Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include "permutation.h"
#include "trsm.h"
#include "ple.h"
#include "newton_john.h"

rci_t mzed_ple_naive(mzed_t *A, mzp_t *P, mzp_t *Q) {
	rci_t col_pos = 0;
	rci_t row_pos = 0;
	word tmp = 0;
	const gf2e *ff = A->finite_field;
	rci_t i,j;
	int found = 0;

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

			if(j+1 < A->ncols) {
				mzed_rescale_row(A, row_pos, j+1, gf2e_inv(ff, tmp));

				for(rci_t l=row_pos+1; l<A->nrows; l++) {
					if ((tmp = mzed_read_elem(A,l,j)))
						mzed_add_multiple_of_row(A, l, A, row_pos, tmp, j+1);
				}
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
	return row_pos;
}

rci_t _mzed_ple(mzed_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff) {
	if (cutoff == 0)
		cutoff = __M4RIE_PLE_CUTOFF;

	if ((A->ncols > m4ri_radix && (gf2e_degree_to_w(A->finite_field) * A->ncols * A->nrows) > cutoff)) {
		mzd_slice_t *a = mzed_slice(NULL, A);
		rci_t r = _mzd_slice_ple(a, P, Q, cutoff);
		mzed_cling(A, a);
		mzd_slice_free(a);
		return r;
	} else {
		return mzed_ple_newton_john(A, P, Q);
	}
}

rci_t _mzd_slice_ple(mzd_slice_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff) {
	const rci_t ncols = A->ncols;
	const rci_t nrows = A->nrows;

	if (cutoff == 0)
		cutoff = __M4RIE_PLE_CUTOFF;

	if (ncols <= m4ri_radix || (gf2e_degree_to_w(A->finite_field) * A->ncols * A->nrows) <= cutoff) {
		mzed_t *Abar = mzed_cling(NULL, A);
		rci_t r = mzed_ple_newton_john(Abar, P, Q);
		mzed_slice(A, Abar);
		mzed_free(Abar);
		return r;
	}

	/*                     n1
	 *   ------------------------------------------
	 *   | A0              | A1                   |
	 *   ------------------------------------------
	 */

	rci_t n1 = (((ncols - 1) / m4ri_radix + 1) >> 1) * m4ri_radix;

	mzd_slice_t *A0 = mzd_slice_init_window(A,  0,  0, nrows,    n1);
	mzd_slice_t *A1 = mzd_slice_init_window(A,  0, n1, nrows, ncols);
	mzp_t *P1 = mzp_init_window(P, 0, nrows);
	mzp_t *Q1 = mzp_init_window(Q, 0, A0->ncols);

	rci_t  r1 = _mzd_slice_ple(A0, P1, Q1, cutoff);

	/*           r1           n1
	 *   ------------------------------------------
	 *   | A00    |           | A01               |
	 * r1------------------------------------------ 
	 *   | A10    |           | A11               |
	 *   ------------------------------------------
	 */

	mzd_slice_t *A00 = mzd_slice_init_window(A,  0,  0, r1, r1);
	mzd_slice_t *A10 = mzd_slice_init_window(A, r1,  0, nrows, r1);
	mzd_slice_t *A01 = mzd_slice_init_window(A,  0, n1, r1, ncols);
	mzd_slice_t *A11 = mzd_slice_init_window(A, r1, n1, nrows, ncols);

	if (r1) {
		/* Computation of the Schur complement */
		mzd_slice_apply_p_left(A1, P1);
		mzd_slice_trsm_lower_left(A00, A01);
		mzd_slice_addmul(A11, A10, A01);
	}
	mzp_free_window(P1);
	mzp_free_window(Q1);

	mzp_t *P2 = mzp_init_window(P, r1, nrows);
	mzp_t *Q2 = mzp_init_window(Q, n1, ncols);

	rci_t r2 = _mzd_slice_ple(A11, P2, Q2, cutoff);

	/*           n
	 *   -------------------
	 *   |      A0b        |
	 *   r1-----------------
	 *   |      A1b        |
	 *   -------------------
	 */

	/* Update A10 */
	mzd_slice_apply_p_left(A10, P2);

	/* Update P */
	for (rci_t i = 0; i < nrows - r1; ++i)
		P2->values[i] += r1;

	/* Update the A0b block (permutation + rotation) */
	for(rci_t i=0, j=n1; j < ncols; ++i, ++j)
		Q2->values[i] += n1;
	for(rci_t i=n1, j = r1; i < n1 + r2; ++i, ++j)
		Q->values[j] = Q->values[i];

	_mzd_slice_compress_l(A, r1, n1, r2);

	mzp_free_window(Q2);
	mzp_free_window(P2);

	mzd_slice_free_window(A0);
	mzd_slice_free_window(A1);
	mzd_slice_free_window(A00);
	mzd_slice_free_window(A01);
	mzd_slice_free_window(A10);
	mzd_slice_free_window(A11);

	return r1 + r2;
}


rci_t _mzd_slice_pluq(mzd_slice_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff) {
	rci_t r = _mzd_slice_ple(A, P, Q, cutoff);
	if(r && r < A->nrows) {
		mzd_slice_t *A0 = mzd_slice_init_window(A, 0, 0, r, A->ncols);
		mzd_slice_apply_p_right_trans_tri(A0, Q);
		mzd_slice_free_window(A0);
	} else {
		mzd_slice_apply_p_right_trans_tri(A, Q);
	}
	return r;
}
