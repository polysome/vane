/**
 * \file m4ri_functions.h
 *
 * \brief Utility functions for mzd_t types
 *
 * \note Some of these functions might be moved M4RI in the future.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_M4RI_FUNCTIONS_H
#define M4RIE_M4RI_FUNCTIONS_H

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

#include <m4ri/m4ri.h>
#include <stdarg.h>

/**
 * mzd_read_bits with assumption that all bits are in the same word
 */

static inline word __mzd_read_bits(const mzd_t *M, const rci_t x, const rci_t y, const rci_t n) {
	int const spot   = y % m4ri_radix;
	wi_t const block = y / m4ri_radix;
	int const spill = spot + n - m4ri_radix;
	word temp = M->rows[x][block] << -spill;
	return temp >> (m4ri_radix - n);
}

/**
 * mzd_xor_bits with assumption that all bits are in the same word
 */

static inline void __mzd_xor_bits(const mzd_t *M, const rci_t x, const rci_t y, const rci_t n, word values) {
	int const spot   = y % m4ri_radix;
	wi_t const block = y / m4ri_radix;
	M->rows[x][block] ^= values << spot;
}

/**
 * mzd_clear_bits with assumption that all bits are in the same word
 */

static inline void __mzd_clear_bits(const mzd_t *M, const rci_t x, const rci_t y, const rci_t n) {
	word values = m4ri_ffff >> (m4ri_radix - n);
	int const spot   = y % m4ri_radix;
	wi_t const block = y / m4ri_radix;
	M->rows[x][block] &= ~(values << spot);
}

/**
 * \brief Add n elements to A
 *
 * A += B[0] + ... + B[n-1]
 *
 * \param A Matrix
 * \param n Number of elements in list
 * \param ... Matrices
 *
 * \ingroup Utility
 */

static inline mzd_t *mzd_sum(mzd_t *A, const int n, ...) {
	assert(n>1);
	va_list b_list;
	va_start( b_list, n );

	mzd_add(A, va_arg(b_list, mzd_t *), va_arg(b_list, mzd_t *));

	for( int i = 0 ; i < n-2; i++ ) {
		mzd_t *B = va_arg(b_list, mzd_t *);
		mzd_add(A, A, B);
	}

	va_end( b_list );
	return A;
}

#endif //M4RIE_M4RI_FUNCTIONS_H
