/**
 * \file mzd_poly.h
 *
 * \brief Matrices over \GF2[x]
 *
 * @warning This code is experimental.
 */

#ifndef M4RIE_MZD_POLY_H
#define M4RIE_MZD_POLY_H

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
#include "mzd_ptr.h"
#include "gf2x.h"
#include "blm.h"

/**
 * \brief will be the data type for matrices over \GF2[x] in the future
 */

typedef struct {
	mzd_t **x;   /**< Coefficients. */
	rci_t nrows; /**< Number of rows. */
	rci_t ncols; /**< Number of columns. */
	deg_t depth;   /**< Degree +1      */
} mzd_poly_t;

/**
 * \brief C += (A+B)*x^offset
 *
 * \param C Target polynomial.
 * \param A Source polynomial.
 * \param B Source polynomial.
 * \param offset The result is shifted offset entries upwards.
 * 
 * \ingroup Addition
 *
 * \warning No bounds checks are performed.
 */

static inline mzd_poly_t *_mzd_poly_add(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B, unsigned int offset) {
	_mzd_ptr_add(C->x+offset, (const mzd_t**)A->x, (const mzd_t**)B->x, A->depth);
	return C;
}

/**
 * \brief C += (A+B)
 *
 * \param C Target polynomial.
 * \param A Source polynomial.
 * \param B Source polynomial.
 * 
 * \ingroup Addition
 */

static inline mzd_poly_t *mzd_poly_add(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B) {
	assert(C->depth >= A->depth && A->depth == B->depth);
	return _mzd_poly_add(C, A, B, 0);
}

/**
 * \brief Create a new polynomial of degree d with m x n matrices as coefficients.
 *
 * \param d Degree.
 * \param m Number of rows.
 * \param n Number of columns.
 *
 * \ingroup Constructions
 */

static inline mzd_poly_t *mzd_poly_init(const deg_t d, const rci_t m, const rci_t n) {
	mzd_poly_t *A = (mzd_poly_t*)m4ri_mm_malloc(sizeof(mzd_poly_t));
	A->x = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*(d+1));

	A->nrows = m;
	A->ncols = n;
	A->depth = d+1;

	for(int i=0; i<A->depth; i++)
		A->x[i] = mzd_init(m,n);
	return A;
}

/**
 * \brief Free polynomial A
 *
 * \param A Polynomial.
 *
 * \ingroup Constructions
 */

static inline void mzd_poly_free(mzd_poly_t *A) {
	for(int i=0; i<A->depth; i++)
		mzd_free(A->x[i]);
	m4ri_mm_free(A->x);
	m4ri_mm_free(A);
}

/**
 * \brief change depth of A to new_depth.
 *
 * \param A Polynomial.
 * \param new_depth New depth (may be <,=,> than current depth).
 *
 * \ingroup Constructions
 */

static inline mzd_poly_t *_mzd_poly_adapt_depth(mzd_poly_t *A, const deg_t new_depth) {
	if (new_depth < A->depth) {
		for(int i=new_depth; i<A->depth; i++) {
			mzd_free(A->x[i]);
			A->x[i] = NULL;
		}
	} else {
		for(int i=A->depth; i<new_depth; i++) {
			A->x[i] = mzd_init(A->nrows,A->ncols);
		}
	}
	A->depth = new_depth;
	return A;
}

/**
 * \brief C += A*B using naive polynomial multiplication
 *
 * \param C Target polynomial.
 * \param A Source polynomial.
 * \param B Source polynomial.
 *
 * \ingroup Multiplication
 */

static inline mzd_poly_t *_mzd_poly_addmul_naive(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B) {
	if (C == NULL)
		C = mzd_poly_init(A->depth+B->depth-1, A->nrows, B->ncols);

	for(unsigned int i=0; i<A->depth; i++) {
		for(unsigned int j=0; j<B->depth; j++) {
			mzd_addmul(C->x[i+j], A->x[i], B->x[j], 0);
		}
	}
	return C;
}

/**
 * \brief C += A*B using Karatsuba multiplication on balanced inputs
 *
 * \param C Target polynomial.
 * \param A Source polynomial.
 * \param B Source polynomial.
 *
 * \ingroup Multiplication
 */


static inline mzd_poly_t *_mzd_poly_addmul_karatsubs_balanced(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B) {
	assert(A->depth == B->depth);

	if (C == NULL)
		C = mzd_poly_init(A->depth+B->depth-1, A->nrows, B->ncols);
	switch(A->depth) {
		case 0:
			m4ri_die("depth 0: seriously?");
		case  1: mzd_addmul(C->x[0], A->x[0], B->x[0], 0); break;
		case  2:  _mzd_ptr_addmul_karatsuba2(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case  3:  _mzd_ptr_addmul_karatsuba3(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case  4:  _mzd_ptr_addmul_karatsuba4(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case  5:  _mzd_ptr_addmul_karatsuba5(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case  6:  _mzd_ptr_addmul_karatsuba6(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case  7:  _mzd_ptr_addmul_karatsuba7(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case  8:  _mzd_ptr_addmul_karatsuba8(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case  9:  _mzd_ptr_addmul_karatsuba9(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case 10: _mzd_ptr_addmul_karatsuba10(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case 11: _mzd_ptr_addmul_karatsuba11(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case 12: _mzd_ptr_addmul_karatsuba12(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case 13: _mzd_ptr_addmul_karatsuba13(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case 14: _mzd_ptr_addmul_karatsuba14(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case 15: _mzd_ptr_addmul_karatsuba15(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		case 16: _mzd_ptr_addmul_karatsuba16(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
		default:
				 m4ri_die("Not implemented\n");
	}
	return C;
}

/**
 * \brief C += A*B by applying the bilinear maps f, i.e. f->H*((f->F*A) x (f->G*B)).
 */

static inline mzd_poly_t *_mzd_poly_addmul_blm(mzd_poly_t *C, mzd_poly_t *A, mzd_poly_t *B, const blm_t *f) {
	assert(f!=NULL);
	assert(f->F->ncols == A->depth && f->G->ncols == B->depth);

	if (C == NULL)
		C = mzd_poly_init(A->depth+B->depth-1, A->nrows, B->ncols);

	_mzd_ptr_apply_blm(C->x, (const mzd_t**)A->x, (const mzd_t**)B->x, f);
	return C;
}


/**
 * \brief C += A*B using the Chinese Remainder Theorem.
 */

static inline mzd_poly_t *_mzd_poly_addmul_crt(mzd_poly_t *C, mzd_poly_t *A, mzd_poly_t *B) {
	int *p = crt_init(A->depth, B->depth);
	blm_t *f = blm_init_crt(NULL, A->depth, B->depth, p, 1);
	_mzd_poly_addmul_blm(C, A, B, f);
	blm_free(f);
	m4ri_mm_free(p);
	return C;
}

/**
 * \brief C += A*B using arithmetic in GF(2^log2(d)) if C has degree d.
 */

mzd_poly_t *_mzd_poly_addmul_ext1(mzd_poly_t *C, mzd_poly_t *A, mzd_poly_t *B);

/**
 * \brief Return -1,0,1 if if A < B, A == B or A > B respectively.
 *
 * \param A Matrix.
 * \param B Matrix.
 *
 * \note This comparison is not well defined (except for !=0) mathematically and relatively
 * arbitrary.
 *
 * \ingroup Comparison
 */

static inline int mzd_poly_cmp(mzd_poly_t *A, mzd_poly_t *B) {
	int r = 0;
	if ((A->depth != B->depth) ) {
		if (A->depth < B->depth)
			return -1;
		else 
			return 1;
	}
	for(int i=0; i<A->depth; i++)
		r |= mzd_cmp(A->x[i],B->x[i]);
	return r;
}

/**
 * \brief Fill matrix A with random elements.
 *
 * \param A Matrix
 *
 * \todo Allow the user to provide a RNG callback.
 *
 * \ingroup Assignment
 */

static inline void mzd_poly_randomize(mzd_poly_t *A) {
	for(int i=0; i<A->depth; i++)
		mzd_randomize(A->x[i]);
}

#endif //M4RIE_MZD_POLY_H
