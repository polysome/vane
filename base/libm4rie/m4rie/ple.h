/**
 * \file ple.h
 * \brief PLE decomposition: \f$ L \cdot E = P \cdot A\f$.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_PLE_H
#define M4RIE_PLE_H

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

#include <m4ri/m4ri.h>
#include <m4rie/mzed.h>
#include <m4rie/mzd_slice.h>
#include <m4rie/conversion.h>

/**
 * \brief PLE decomposition: \f$ L \cdot E = P \cdot A \f$.
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- an echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function uses naive cubic PLE decomposition depending on the
 * size of the underlying field.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 *
 * \ingroup PLE
 *
 * \sa mzed_ple_newton_john() mzed_ple()
 */

rci_t mzed_ple_naive(mzed_t *A, mzp_t *P, mzp_t *Q);

/**
 * \brief PLE decomposition: \f$ L \cdot E = P \cdot A \f$.
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- an echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function uses either asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication or naive cubic PLE
 * decomposition depending on the size of the underlying field. If
 * asymptotically fast PLE decomposition is used, then the algorithm
 * switches to mzed_ple_newton_john if e * ncols * nrows is <= cutoff
 * where e is the exponent of the finite field.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 * \param cutoff Integer
 *
 * \ingroup PLE
 *
 * \sa mzed_ple_naive() mzed_ple_newton_john() mzed_ple()
 */

rci_t _mzd_slice_ple(mzd_slice_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff);

/**
 * \brief PLE decomposition: \f$ L \cdot E = P \cdot A \f$.
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- an echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function implements asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 *
 * \ingroup PLE
 *
 * \sa mzed_ple_naive() mzed_ple_newton_john() _mzd_slice_ple()
 */

static inline rci_t mzd_slice_ple(mzd_slice_t *A, mzp_t *P, mzp_t *Q) {
	assert(P->length == A->nrows);
	assert(Q->length == A->ncols);
	return _mzd_slice_ple(A, P, Q, 0);
}

/**
 * \brief PLUQ decomposition: \f$ L \cdot U \cdot Q = P \cdot A\f$.
 *
 * This function implements asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication. From PLE the PLUQ
 * decomposition is then obtained.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 * \param cutoff Crossover to base case if mzed_t::w * mzed_t::ncols * mzed_t::nrows < cutoff.
 *
 * \ingroup PLE
 */

rci_t _mzd_slice_pluq(mzd_slice_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff);

/**
 * \brief PLUQ decomposition: \f$ L \cdot U \cdot Q = P \cdot A\f$.
 *
 * This function implements asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication. From PLE the PLUQ
 * decomposition is then obtained.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 *
 * \ingroup PLE
 */

static inline rci_t mzd_slice_pluq(mzd_slice_t *A, mzp_t *P, mzp_t *Q) {
	assert(P->length == A->nrows);
	assert(Q->length == A->ncols);
	return _mzd_slice_pluq(A, P, Q, 0);
}


/**
 * \brief PLE decomposition: \f$ L \cdot E = P \cdot A \f$.
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- an echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function uses either asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication or naive cubic PLE
 * decomposition depending on the size of the underlying field. If
 * asymptotically fast PLE decomposition is used, then the algorithm
 * switches to mzed_ple_newton_john if e * ncols * nrows is <= cutoff
 * where e is the exponent of the finite field.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 * \param cutoff Integer >= 0
 *
 * \ingroup PLE
 *
 * \sa mzed_ple_naive() mzed_ple_newton_john() _mzed_ple()
 */

rci_t _mzed_ple(mzed_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff);

/**
 * Default crossover to PLE base case (Newton-John based).
 */

#define __M4RIE_PLE_CUTOFF (__M4RI_CPU_L2_CACHE<<2)

/**
 * \brief PLE decomposition: \f$ L \cdot E = P \cdot A \f$.
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- an echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function uses either asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication or naive cubic PLE
 * decomposition depending on the size of the underlying field.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 *
 * \ingroup PLE
 *
 * \sa mzed_ple_naive() mzed_ple_newton_john() _mzed_ple()
 *
 */

static inline rci_t mzed_ple(mzed_t *A, mzp_t *P, mzp_t *Q) {
	return _mzed_ple(A, P, Q, __M4RIE_PLE_CUTOFF);
}

#endif //M4RIE_PLE_H
