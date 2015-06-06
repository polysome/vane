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

#include "config.h"

#include <m4ri/misc.h>
#include <m4ri/mzd.h>
#include <m4ri/m4ri_config.h>

#include "mzed.h"
#include "newton_john.h"
#include "mzd_slice.h"
#include "strassen.h"

#define CLOSER(a,b,target) (abs((long)a-(long)target)<abs((long)b-(long)target))

mzed_t *mzed_mul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff) {
	C = _mzed_mul_init(C, A, B, TRUE);
	return _mzed_mul_strassen(C, A, B, cutoff);
}

mzed_t *mzed_addmul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff) {
	C = _mzed_mul_init(C, A, B, FALSE);
	return _mzed_addmul_strassen(C, A, B, cutoff);
}

mzed_t *_mzed_mul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff) {
	if((C->nrows| C->ncols) == 0)
		return C;

	rci_t m = A->nrows;
	rci_t k = A->ncols;
	rci_t n = B->ncols;

	/* handle case first, where the input matrices are too small already */
	if (CLOSER(m, m/2, cutoff) || CLOSER(k, k/2, cutoff) || CLOSER(n, n/2, cutoff)) {
		/* we copy the matrix first since it is only constant memory
		   overhead and improves data locality, if you remove it make sure
		   there are no speed regressions */
		/* C = _mzd_mul_m4rm(C, A, B, 0, TRUE); */
		mzed_t *Cbar = mzed_init(C->finite_field, C->nrows, C->ncols);
		_mzed_mul_newton_john(Cbar, A, B);

		mzed_copy(C, Cbar);
		mzed_free(Cbar);
		return C;
	}

	rci_t mmm = m/2;
	rci_t kkk = k/2;
	rci_t nnn = n/2;

	mmm = (mmm - mmm%(m4ri_radix/A->w));
	kkk = (kkk - kkk%(m4ri_radix/A->w));
	nnn = (nnn - nnn%(m4ri_radix/A->w));

	/*         |A |   |B |   |C |
	 * Compute |  | x |  | = |  | */

	mzed_t *A11 = mzed_init_window(A,   0,   0,   mmm,   kkk);
	mzed_t *A12 = mzed_init_window(A,   0, kkk,   mmm, 2*kkk);
	mzed_t *A21 = mzed_init_window(A, mmm,   0, 2*mmm,   kkk);
	mzed_t *A22 = mzed_init_window(A, mmm, kkk, 2*mmm, 2*kkk);

	mzed_t *B11 = mzed_init_window(B,   0,   0,   kkk,   nnn);
	mzed_t *B12 = mzed_init_window(B,   0, nnn,   kkk, 2*nnn);
	mzed_t *B21 = mzed_init_window(B, kkk,   0, 2*kkk,   nnn);
	mzed_t *B22 = mzed_init_window(B, kkk, nnn, 2*kkk, 2*nnn);

	mzed_t *C11 = mzed_init_window(C,   0,   0,   mmm,   nnn);
	mzed_t *C12 = mzed_init_window(C,   0, nnn,   mmm, 2*nnn);
	mzed_t *C21 = mzed_init_window(C, mmm,   0, 2*mmm,   nnn);
	mzed_t *C22 = mzed_init_window(C, mmm, nnn, 2*mmm, 2*nnn);

	/**
	 * \note See Marco Bodrato; "A Strassen-like Matrix Multiplication
	 * Suited for Squaring and Highest Power Computation";
	 * http://bodrato.it/papres/#CIVV2008 for reference on the used
	 * sequence of operations.
	 */

	/* change this to mzd_init(mmm, MAX(nnn,kkk)) to fix the todo below */
	mzed_t *Wmk = mzed_init(A->finite_field, mmm, kkk);
	mzed_t *Wkn = mzed_init(A->finite_field, kkk, nnn);

	_mzed_add(Wkn, B22, B12);		 /* Wkn = B22 + B12 */
	_mzed_add(Wmk, A22, A12);		 /* Wmk = A22 + A12 */
	_mzed_mul_strassen(C21, Wmk, Wkn, cutoff);    /* C21 = Wmk * Wkn */

	_mzed_add(Wmk, A22, A21);		 /* Wmk = A22 - A21 */
	_mzed_add(Wkn, B22, B21);		 /* Wkn = B22 - B21 */
	_mzed_mul_strassen(C22, Wmk, Wkn, cutoff);    /* C22 = Wmk * Wkn */

	_mzed_add(Wkn, Wkn, B12);		   /* Wkn = Wkn + B12 */
	_mzed_add(Wmk, Wmk, A12);		   /* Wmk = Wmk + A12 */
	_mzed_mul_strassen(C11, Wmk, Wkn, cutoff); /* C11 = Wmk * Wkn */

	_mzed_add(Wmk, Wmk, A11);		 /* Wmk = Wmk - A11 */
	_mzed_mul_strassen(C12, Wmk, B12, cutoff);    /* C12 = Wmk * B12 */
	_mzed_add(C12, C12, C22);		 /* C12 = C12 + C22 */

	mzed_free(Wmk);
	Wmk = mzed_mul_strassen(NULL, A12, B21, cutoff);/*Wmk = A12 * B21 */

	_mzed_add(C11, C11, Wmk);		  /* C11 = C11 + Wmk */
	_mzed_add(C12, C11, C12);		  /* C12 = C11 - C12 */
	_mzed_add(C11, C21, C11);		  /* C11 = C21 - C11 */
	_mzed_add(Wkn, Wkn, B11);		  /* Wkn = Wkn - B11 */
	_mzed_mul_strassen(C21, A21, Wkn, cutoff);     /* C21 = A21 * Wkn */
	mzed_free(Wkn);

	_mzed_add(C21, C11, C21);		  /* C21 = C11 - C21 */
	_mzed_add(C22, C22, C11);		  /* C22 = C22 + C11 */
	_mzed_mul_strassen(C11, A11, B11, cutoff);     /* C11 = A11 * B11 */

	_mzed_add(C11, C11, Wmk);		  /* C11 = C11 + Wmk */

	/* clean up */
	mzed_free_window(A11); mzed_free_window(A12);
	mzed_free_window(A21); mzed_free_window(A22);

	mzed_free_window(B11); mzed_free_window(B12);
	mzed_free_window(B21); mzed_free_window(B22);

	mzed_free_window(C11); mzed_free_window(C12);
	mzed_free_window(C21); mzed_free_window(C22);

	mzed_free(Wmk);

	/* deal with rest */
	nnn*=2;
	if (n > nnn) {
		/*         |AA|   | B|   | C|
		 * Compute |AA| x | B| = | C| */
		mzed_t *B_last_col = mzed_init_window(B, 0, nnn, k, n);
		mzed_t *C_last_col = mzed_init_window(C, 0, nnn, m, n);
		mzed_set_ui(C_last_col, 0);
		_mzed_mul_newton_john(C_last_col, A, B_last_col);
		mzed_free_window(B_last_col);
		mzed_free_window(C_last_col);
	}
	mmm*=2;
	if (m > mmm) {
		/*         |  |   |B |   |  |
		 * Compute |AA| x |B | = |C | */
		mzed_t *A_last_row = mzed_init_window(A, mmm, 0, m, k);
		mzed_t *B_first_col= mzed_init_window(B,   0, 0, k, nnn);
		mzed_t *C_last_row = mzed_init_window(C, mmm, 0, m, nnn);
		mzed_set_ui(C_last_row, 0);
		_mzed_mul_newton_john(C_last_row, A_last_row, B_first_col);
		mzed_free_window(A_last_row);
		mzed_free_window(B_first_col);
		mzed_free_window(C_last_row);
	}
	kkk*=2;
	if (k > kkk) {
		/* Add to  |  |   | B|   |C |
		 * result  |A | x |  | = |  | */
		mzed_t *A_last_col = mzed_init_window(A,   0, kkk, mmm, k);
		mzed_t *B_last_row = mzed_init_window(B, kkk,   0,   k, nnn);
		mzed_t *C_bulk = mzed_init_window(C, 0, 0, mmm, nnn);
		_mzed_mul_newton_john(C_bulk, A_last_col, B_last_row);
		mzed_free_window(A_last_col);
		mzed_free_window(B_last_row);
		mzed_free_window(C_bulk);
	}

	return C;
}

mzed_t *_mzed_addmul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff) {
	if((C->nrows| C->ncols) == 0)
		return C;

	rci_t m = A->nrows;
	rci_t k = A->ncols;
	rci_t n = B->ncols;

	/* handle case first, where the input matrices are too small already */
	if (CLOSER(m, m/2, cutoff) || CLOSER(k, k/2, cutoff) || CLOSER(n, n/2, cutoff)) {
		/* we copy the matrix first since it is only constant memory
		   overhead and improves data locality, if you remove it make sure
		   there are no speed regressions */
		/* C = _mzd_mul_m4rm(C, A, B, 0, TRUE); */
		mzed_t *Cbar = mzed_copy(NULL, C);
		_mzed_mul_newton_john(Cbar, A, B);

		mzed_copy(C, Cbar);
		mzed_free(Cbar);
		return C;
	}

	rci_t mmm = m/2;
	rci_t kkk = k/2;
	rci_t nnn = n/2;

	mmm = (mmm - mmm%(m4ri_radix/A->w));
	kkk = (kkk - kkk%(m4ri_radix/A->w));
	nnn = (nnn - nnn%(m4ri_radix/A->w));

	/*         |A |   |B |   |C |
	 * Compute |  | x |  | = |  | */

	mzed_t *A11 = mzed_init_window(A,   0,   0,   mmm,   kkk);
	mzed_t *A12 = mzed_init_window(A,   0, kkk,   mmm, 2*kkk);
	mzed_t *A21 = mzed_init_window(A, mmm,   0, 2*mmm,   kkk);
	mzed_t *A22 = mzed_init_window(A, mmm, kkk, 2*mmm, 2*kkk);

	mzed_t *B11 = mzed_init_window(B,   0,   0,   kkk,   nnn);
	mzed_t *B12 = mzed_init_window(B,   0, nnn,   kkk, 2*nnn);
	mzed_t *B21 = mzed_init_window(B, kkk,   0, 2*kkk,   nnn);
	mzed_t *B22 = mzed_init_window(B, kkk, nnn, 2*kkk, 2*nnn);

	mzed_t *C11 = mzed_init_window(C,   0,   0,   mmm,   nnn);
	mzed_t *C12 = mzed_init_window(C,   0, nnn,   mmm, 2*nnn);
	mzed_t *C21 = mzed_init_window(C, mmm,   0, 2*mmm,   nnn);
	mzed_t *C22 = mzed_init_window(C, mmm, nnn, 2*mmm, 2*nnn);

	/**
	 * \note See Marco Bodrato; "A Strassen-like Matrix Multiplication
	 * Suited for Squaring and Highest Power Computation";
	 * http://bodrato.it/papres/#CIVV2008 for reference on the used
	 * sequence of operations.
	 */

	mzed_t *S = mzed_init(A->finite_field, mmm, kkk);
	mzed_t *T = mzed_init(A->finite_field, kkk, nnn);
	mzed_t *U = mzed_init(A->finite_field, mmm, nnn);

	_mzed_add(S, A22, A21);                   /* 1  S = A22 - A21       */
	_mzed_add(T, B22, B21);                   /* 2  T = B22 - B21       */
	_mzed_mul_strassen(U, S, T, cutoff); /* 3  U = S*T             */
	_mzed_add(C22, U, C22);                   /* 4  C22 = U + C22       */
	_mzed_add(C12, U, C12);                   /* 5  C12 = U + C12       */

	_mzed_mul_strassen(U, A12, B21, cutoff); /* 8  U = A12*B21         */
	_mzed_add(C11, U, C11);                       /* 9  C11 = U + C11       */

	_mzed_addmul_strassen(C11, A11, B11, cutoff); /* 11 C11 = A11*B11 + C11 */

	_mzed_add(S, S, A12);                     /* 6  S = S - A12         */
	_mzed_add(T, T, B12);                     /* 7  T = T - B12         */
	_mzed_addmul_strassen(U, S, T, cutoff); /* 10 U = S*T + U         */
	_mzed_add(C12, C12, U);                   /* 15 C12 = U + C12       */

	_mzed_add(S, A11, S);                     /* 12 S = A11 - S         */
	_mzed_addmul_strassen(C12, S, B12, cutoff);   /* 14 C12 = S*B12 + C12   */

	_mzed_add(T, B11, T);                     /* 13 T = B11 - T         */
	_mzed_addmul_strassen(C21, A21, T, cutoff);   /* 16 C21 = A21*T + C21   */

	_mzed_add(S, A22, A12);                   /* 17 S = A22 + A21       */
	_mzed_add(T, B22, B12);                   /* 18 T = B22 + B21       */
	_mzed_addmul_strassen(U, S, T, cutoff);       /* 19 U = U - S*T         */
	_mzed_add(C21, C21, U);                   /* 20 C21 = C21 - U3      */
	_mzed_add(C22, C22, U);                   /* 21 C22 = C22 - U3      */

	/* clean up */
	mzed_free_window(A11); mzed_free_window(A12);
	mzed_free_window(A21); mzed_free_window(A22);

	mzed_free_window(B11); mzed_free_window(B12);
	mzed_free_window(B21); mzed_free_window(B22);

	mzed_free_window(C11); mzed_free_window(C12);
	mzed_free_window(C21); mzed_free_window(C22);

	mzed_free(S);
	mzed_free(T);
	mzed_free(U);


	/* deal with rest */
	nnn*=2;
	if (n > nnn) {
		/*         |AA|   | B|   | C|
		 * Compute |AA| x | B| = | C| */
		mzed_t const *B_last_col = mzed_init_window(B, 0, nnn, k, n); 
		mzed_t *C_last_col = mzed_init_window(C, 0, nnn, m, n);
		_mzed_mul_newton_john(C_last_col, A, B_last_col);
		mzed_free_window((mzed_t*)B_last_col);
		mzed_free_window((mzed_t*)C_last_col);
	}
	mmm*=2;
	if (m > mmm) {
		/*         |  |   |B |   |  |
		 * Compute |AA| x |B | = |C | */
		mzed_t const *A_last_row = mzed_init_window(A, mmm, 0, m, k);
		mzed_t const *B_first_col= mzed_init_window(B,   0, 0, k, nnn);
		mzed_t *C_last_row = mzed_init_window(C, mmm, 0, m, nnn);
		_mzed_mul_newton_john(C_last_row, A_last_row, B_first_col);
		mzed_free_window((mzed_t*)A_last_row);
		mzed_free_window((mzed_t*)B_first_col);
		mzed_free_window(C_last_row);
	}
	kkk*=2;
	if (k > kkk) {
		/* Add to  |  |   | B|   |C |
		 * result  |A | x |  | = |  | */
		mzed_t const *A_last_col = mzed_init_window(A,   0, kkk, mmm, k);
		mzed_t const *B_last_row = mzed_init_window(B, kkk,   0,   k, nnn);
		mzed_t *C_bulk = mzed_init_window(C, 0, 0, mmm, nnn);
		_mzed_mul_newton_john(C_bulk, A_last_col, B_last_row);
		mzed_free_window((mzed_t*)A_last_col);
		mzed_free_window((mzed_t*)B_last_row);
		mzed_free_window(C_bulk);
	}

	return C;
}

rci_t _mzed_strassen_cutoff(const mzed_t *C, const mzed_t *A, const mzed_t *B) {
	rci_t cutoff;

	/* it seems most of it is cache bound: 2 matrix * (n^2 *w / 8 ) <= L2  */

	switch(A->finite_field->degree) {

		case 2:
			cutoff = MIN(((int)sqrt((double)(4*__M4RI_CPU_L2_CACHE)))/2,4096);
			break;
		case  3:
		case  4:
		case  5:
		case  6:
		case  7:
		case  8:
			cutoff = MIN(((int)sqrt((double)(4*__M4RI_CPU_L2_CACHE/A->w))),4096);
			break;

		case  9:
			/* on redhawk 2048 is much better, sage.math 1204 wins **/
			cutoff = 2048;
			break;

		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
			cutoff = 4096;
			break;

		default:
			cutoff = 1024;
			break;
	}

	if (cutoff < 2*__M4RI_TWOPOW(A->finite_field->degree))
		cutoff = 2*__M4RI_TWOPOW(A->finite_field->degree);
	return cutoff;
}
