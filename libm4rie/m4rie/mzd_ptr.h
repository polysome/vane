#ifndef M4RIE_MZD_PTR_H
#define M4RIE_MZD_PTR_H

/******************************************************************************
 *
 *            M4RIE: Linear Algebra over GF(2^e)
 *
 *    Copyright (C) 2013 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include <stdarg.h>
#include <m4rie/gf2e.h>
#include <m4ri/djb.h>

/**
 * \brief Add A to coefficient of X^t but perform modular reductions on the fly.
 *
 * (A + X^t % minpoly)
 *
 * \param ff Finite field (must not be NULL)
 * \param A Matrix
 * \param X Matrix list
 * \param t Integer >= 0 (degree)
 */

static inline void _mzd_ptr_add_modred(const gf2e *ff, const mzd_t *A, mzd_t **X, const int t) {
	if (mzd_is_zero(A))
		return;

	if ((ff==NULL) || (t < ff->degree)) {
		mzd_add(X[t], X[t], A);
		return;
	}

	word pow_gen = ff->pow_gen[t];

	for(int i=0; i<ff->degree; i++) {
		if (pow_gen & (1<<i))
			mzd_add(X[i],X[i],A);
	}
}

/**
 * \brief Add A to n coefficients but perform modular reductions on the fly (if ff != NULL).
 *
 * \param ff Finite field (may be NULL)
 * \param A Matrix
 * \param X Matrix list
 * \param n Integer > 0
 */

static inline mzd_t *_mzd_ptr_add_to_all(const gf2e *ff, mzd_t *A, mzd_t **X, const int n, ...) {
	va_list b_list;
	va_start( b_list, n );

	if (ff) 
		for( int i = 0 ; i < n; i++ ) {
			int t = va_arg(b_list, int);
			_mzd_ptr_add_modred(ff, A, X, t);
		}
	else
		for( int i = 0 ; i < n; i++ ) {
			int t = va_arg(b_list, int);
			mzd_add(X[t], X[t], A);
		}

	va_end( b_list );
	return A;
}

static inline void _mzd_ptr_add(mzd_t **c, const mzd_t **a, const mzd_t **b, const deg_t length) {
	switch(length) {
		case 32: mzd_add(c[31], a[31], b[31]);
		case 31: mzd_add(c[30], a[30], b[30]);
		case 30: mzd_add(c[29], a[29], b[29]);
		case 29: mzd_add(c[28], a[28], b[28]);
		case 28: mzd_add(c[27], a[27], b[27]);
		case 27: mzd_add(c[26], a[26], b[26]);
		case 26: mzd_add(c[25], a[25], b[25]);
		case 25: mzd_add(c[24], a[24], b[24]);
		case 24: mzd_add(c[23], a[23], b[23]);
		case 23: mzd_add(c[22], a[22], b[22]);
		case 22: mzd_add(c[21], a[21], b[21]);
		case 21: mzd_add(c[20], a[20], b[20]);
		case 20: mzd_add(c[19], a[19], b[19]);
		case 19: mzd_add(c[18], a[18], b[18]);
		case 18: mzd_add(c[17], a[17], b[17]);
		case 17: mzd_add(c[16], a[16], b[16]);
		case 16: mzd_add(c[15], a[15], b[15]);
		case 15: mzd_add(c[14], a[14], b[14]);
		case 14: mzd_add(c[13], a[13], b[13]);
		case 13: mzd_add(c[12], a[12], b[12]);
		case 12: mzd_add(c[11], a[11], b[11]);
		case 11: mzd_add(c[10], a[10], b[10]);
		case 10: mzd_add(c[ 9], a[ 9], b[ 9]);
		case  9: mzd_add(c[ 8], a[ 8], b[ 8]);
		case  8: mzd_add(c[ 7], a[ 7], b[ 7]);
		case  7: mzd_add(c[ 6], a[ 6], b[ 6]);
		case  6: mzd_add(c[ 5], a[ 5], b[ 5]);
		case  5: mzd_add(c[ 4], a[ 4], b[ 4]);
		case  4: mzd_add(c[ 3], a[ 3], b[ 3]);
		case  3: mzd_add(c[ 2], a[ 2], b[ 2]);
		case  2: mzd_add(c[ 1], a[ 1], b[ 1]);
		case  1: mzd_add(c[ 0], a[ 0], b[ 0]);
		case  0:
				 break;
		default:
				 for(deg_t i=0; i<length; i++)
					 mzd_add(c[ i], a[ i], b[ i]);
	}
}

/** karatsuba.c **/

/**
 * \brief \f$ X += A \cdot B \f$ over \GF4 using 3 multiplications over \GF2 and 2 temporary \GF2 matrices.
 *
 * If no finite field is given, polynomial arithmetic with polynomials of degree 1 is performed. In
 * this case, X is expected to have at least length 3. If a finite field is given, then C is
 * expected to have at least length 2.
 *
 * The formula was taken from Peter L. Montgomery. "Five, Six, and
 * Seven-Term Karatsuba-Like Formulae" in IEEE TRANSACTIONS ON
 * COMPUTERS, VOL. 54, NO. 3, MARCH 2005/
 * 
 * \param ff Finite Field, may be NULL for polynomial arithmetic.
 * \param X Preallocated return matrix, of length >= 2 (ff != NULL) or >=3 (ff == NULL)
 * \param A Input matrix A, preallocated of length >= 2.
 * \param B Input matrix B, preallocated of length >= 2.
 *
 * \sa _mzd_ptr_addmul_karatsuba()
 *
 * \ingroup Multiplication
 */

void _mzd_ptr_addmul_karatsuba2(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF8 using 6 multiplications over \GF2 and 3 temporary \GF2 matrices..
 */

void _mzd_ptr_addmul_karatsuba3(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF16 using 9 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba4(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ C += A \cdot B \f$ over \GF32 using 13 multiplications over \GF2 and 3 temporary \GF2 matrices..
 */

void _mzd_ptr_addmul_karatsuba5(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ C += A \cdot B \f$ over \GF64 using 17 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba6(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF128 using 22 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba7(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF256 using 27 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba8(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF512 using 31 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba9(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF1024 using 36 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba10(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF2048 using 40 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba11(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF4096 using 45 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba12(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF8192 using 49 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba13(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF16384 using 55 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba14(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF32768 using 60 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba15(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);

/**
 * \brief \f$ X += A \cdot B \f$ over \GF65536 using 64 multiplications over \GF2 and 3 temporary \GF2 matrices.
 */

void _mzd_ptr_addmul_karatsuba16(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B);


/**
 * W = m*V
 *
 * Apply the linear map m to V (considered as a vector) and write the result to W.
 *
 * \param m Linear map
 * \param W Output vector of matrices mod 2.
 * \oaram V Input vector of matrices mod 2.
 */

void djb_apply_mzd_ptr(djb_t *m, mzd_t **W, const mzd_t **V);

#endif //M4RIE_MZD_PTR_H
