/**
 * \file trsm.h
 * \brief Triangular System Solving with Matrices (TRSM).
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef TRSM_H
#define TRSM_H

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

#include <m4rie/mzed.h>
#include <m4rie/mzd_slice.h>

#define MZED_TRSM_CUTOFF 512 /**< Crossover dimension to TRSM base cases */

/**
 * \brief \f$B = U^{-1} \cdot B\f$
 *
 * \param U Upper-triangular matrix (other entries are ignored).
 * \param B Matrix.
 * \param cutoff Crossover dimension to base case.
 *
 * \ingroup Triangular
 */

void _mzed_trsm_upper_left(const mzed_t *U, mzed_t *B, const rci_t cutoff);

/**
 * \brief \f$B = U^{-1} \cdot B\f$
 *
 * \param U Upper-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

void mzed_trsm_upper_left_naive(const mzed_t *U, mzed_t *B);

/**
 * \brief \f$B = U^{-1} \cdot B\f$
 *
 * \param U Upper-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

static inline void mzed_trsm_upper_left(const mzed_t *U, mzed_t *B) {
	_mzed_trsm_upper_left(U, B, MZED_TRSM_CUTOFF);
}

/**
 * \brief \f$B = U^{-1} \cdot B\f$
 *
 * \param U Upper-triangular matrix (other entries are ignored).
 * \param B Matrix.
 * \param cutoff Crossover dimension to base case.
 *
 * \ingroup Triangular
 */

void _mzd_slice_trsm_upper_left(const mzd_slice_t *U, mzd_slice_t *B, const rci_t cutoff);

/**
 * \brief \f$B = U^{-1} \cdot B\f$
 *
 * \param U Upper-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

void mzd_slice_trsm_upper_left_naive(const mzd_slice_t *U, mzd_slice_t *B);

/**
 * \brief \f$B = U^{-1} \cdot B\f$
 *
 * \param U Upper-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

static inline void mzd_slice_trsm_upper_left(const mzd_slice_t *U, mzd_slice_t *B) {
	_mzd_slice_trsm_upper_left(U, B, MZED_TRSM_CUTOFF);
}

/**
 * \brief \f$B = L^{-1} \cdot B\f$
 *
 * \param L Lower-triangular matrix (other entries are ignored).
 * \param B Matrix.
 * \param cutoff Crossover dimension to base case.
 *
 * \ingroup Triangular
 */

void _mzed_trsm_lower_left(const mzed_t *L, mzed_t *B, const rci_t cutoff);

/**
 * \brief \f$B = L^{-1} \cdot B\f$
 *
 * \param L Lower-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

void mzed_trsm_lower_left_naive(const mzed_t *L, mzed_t *B);

/**
 * \brief \f$B = L^{-1} \cdot B\f$
 *
 * \param L Lower-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

static inline void mzed_trsm_lower_left(const mzed_t *L, mzed_t *B) {
	_mzed_trsm_lower_left(L, B, MZED_TRSM_CUTOFF);
}

/**
 * \brief \f$B = L^{-1} \cdot B\f$
 *
 * \param L Lower-triangular matrix (other entries are ignored).
 * \param B Matrix.
 * \param cutoff Crossover dimension to base case.
 *
 * \ingroup Triangular
 */

void _mzd_slice_trsm_lower_left(const mzd_slice_t *L, mzd_slice_t *B, const rci_t cutoff);

/**
 * \brief \f$B = L^{-1} \cdot B\f$
 *
 * \param L Lower-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

void mzd_slice_trsm_lower_left_naive(const mzd_slice_t *L, mzd_slice_t *B);

/**
 * \brief \f$B = L^{-1} \cdot B\f$
 *
 * \param L Lower-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

static inline void mzd_slice_trsm_lower_left(const mzd_slice_t *L, mzd_slice_t *B) {
	_mzd_slice_trsm_lower_left(L, B, MZED_TRSM_CUTOFF);
}




#endif //TRSM_H
