/**
 * \file permutation.h
 * \brief Permutation matrices.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_PERMUTATION_H
#define M4RIE_PERMUTATION_H

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

#include <m4ri/mzp.h>
#include <m4rie/mzed.h>
#include <m4rie/mzd_slice.h>

/**
 * Apply the permutation P to A from the left.
 *
 * This is equivalent to row swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzed_apply_p_left(mzed_t *A, mzp_t const *P) {
	mzd_apply_p_left(A->x, P);
}

/**
 * Apply the permutation P to A from the left but transpose P before.
 *
 * This is equivalent to row swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzed_apply_p_left_trans(mzed_t *A, mzp_t const *P) {
	mzd_apply_p_left_trans(A->x, P);
}

/**
 * Apply the permutation P to A from the right.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzed_apply_p_right(mzed_t *A, mzp_t const *P) {
	if(A->nrows == 0)
		return;
	rci_t const length = MIN(P->length, A->ncols);
	for (rci_t i = length-1; i >= 0; --i) {
		mzed_col_swap(A, i, P->values[i]);
	}
}

/**
 * Apply the permutation P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzed_apply_p_right_trans(mzed_t *A, mzp_t const *P) {
	if(A->nrows == 0)
		return;
	rci_t const length = MIN(P->length, A->ncols);
	for (rci_t i = 0; i < length; ++i) {
		mzed_col_swap(A, i, P->values[i]);
	}
}

/**
 * Apply the permutation P to A from the left.
 *
 * This is equivalent to row swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzd_slice_apply_p_left(mzd_slice_t *A, mzp_t const *P) {
	for(int i=0; i<A->depth; i++) {
		mzd_apply_p_left(A->x[i], P);
	}
}

/**
 * Apply the permutation P to A from the left but transpose P before.
 *
 * This is equivalent to row swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzd_slice_apply_p_left_trans(mzd_slice_t *A, mzp_t const *P) {
	for(int i=0; i<A->depth; i++) {
		mzd_apply_p_left_trans(A->x[i], P);
	}
}

/**
 * Apply the permutation P to A from the right.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzd_slice_apply_p_right(mzd_slice_t *A, mzp_t const *P) {
	for(int i=0; i<A->depth; i++) {
		mzd_apply_p_right(A->x[i], P);
	}
}

/**
 * Apply the permutation P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzd_slice_apply_p_right_trans(mzd_slice_t *A, mzp_t const *P) {
	for(int i=0; i<A->depth; i++) {
		mzd_apply_p_right_trans(A->x[i], P);
	}
}

/**
 * Apply the permutation P to A from the right, but only on the upper
 * the matrix A above the main diagonal.
 *
 * This is equivalent to column swaps walking from 0 to length-1 and
 * is used to compress PLE to PLUQ.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

static inline void mzd_slice_apply_p_right_trans_tri(mzd_slice_t *A, mzp_t const *P) {
	for(int i=0; i<A->depth; i++) {
		mzd_apply_p_right_trans_tri(A->x[i], P);
	}
}


#endif // M4RIE_PERMUTATION_H
