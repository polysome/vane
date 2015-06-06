/**
 * \file blm.h
 *
 * \brief Bilinear Maps on Matrices over GF(2).
 *
 * This is used to realise mzd_poly_t multiplication.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_BLM_H
#define M4RIE_BLM_H

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

#include <m4ri/m4ri.h>
#include "m4rie/gf2e.h"

/**
 * \brief We consider at most polynomials of degree M4RIE_MAX_DEGREE in CRT.
 */

#define M4RIE_CRT_LEN (M4RIE_MAX_DEGREE + 1)

/**
 * \brief Bilinear Maps on Matrices over GF(2).
 *
 * Encodes the bilinear map H*((F*A) x (G*B)) where A,B are vectors of mzd_t, "*" is matrix-vector
 * multiplication and "x" is pointwise multiplication.
 *
 * If a DJB map is not NULL, it will be used instead its matrix representant.
 */

typedef struct {
	mzd_t *H; /*!< final linear map H*/
	djb_t *h; /*!< final linear map H (DJB encoding)*/

	mzd_t *F; /*!< lineatr map on A */
	djb_t *f; /*!< lineatr map on N (DJB encoding) */

	mzd_t *G; /*!< lineatr map on B */
	djb_t *g; /*!< lineatr map on B (DJB encoding) */
} blm_t;

/**
 * costs[i] = cost of multiplying two polynomials of length i over \GF2.
 */

extern const int costs[17];

/**
 * Return the multiplication cost of the multiplication scheme p
 */

static inline int blm_cost_crt(const int p[M4RIE_CRT_LEN]) {
	int cost = costs[p[0]];
	for(deg_t d=1; d<M4RIE_CRT_LEN; d++)
		cost += costs[d] * p[d];
	return cost;
}

/**
 * Find a list of co-prime polynomials p_i such that deg(prod(p_i)) >= f_len*g_len-1.
 *
 * We store the number of polynomials of degree d in p[d]. We store the degree w of (x-infinity)^w
 *  in p[0].
 */

int *crt_init(const deg_t f_len, const deg_t g_len);

/**
 * Compute H, F, G such that vec(c) = H*(F*vec(a) x G*vec(b))
 * with poly(c) = poly(a)*poly(b), 
 * deg(poly(a)) = a_ncols -1, deg(poly(b)) = b_ncols -1 and "x" being pointwise multiplication
 *
 * This is realised by a multi-modular computation modulo the primes up to degree deg (which may not
 * be irreducible polynomials, but merely co-prime).
 */

blm_t *blm_init_crt(const gf2e *ff, const deg_t f_ncols, const deg_t g_ncols, const int *p, int djb);

/**
 * Given F and G compute H.
 *
 * \param ff finite field for modular reduction
 * \param f Bilinear Map with F and G already computed.
 */

blm_t *_blm_finish_polymult(const gf2e *ff, blm_t *f);

/**
 * Free bilinear map f.
 */

void blm_free(blm_t *f);

/**
 * \brief Compile DJB map for f
 *  
 * \param f Bilinear map
 */

blm_t *_blm_djb_compile(blm_t *f);

/**
 * \brief Apply f (stored as a matrix) on A and B, writing to X
 *  
 * \param X Array of matrices
 * \param A Array of matrices 
 * \param B Array of matrices 
 * \param f Bilinear map
 */

void _mzd_ptr_apply_blm_mzd(mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f);

/**
 * \brief Apply f (stored as a DJB map) on A and B, writing to X
 *  
 * \param X Array of matrices
 * \param A Array of matrices 
 * \param B Array of matrices 
 * \param f Bilinear map
 */

void _mzd_ptr_apply_blm_djb(mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f);

/**
 * \brief Apply f on A and B, writing to X
 *  
 * \param X Array of matrices
 * \param A Array of matrices 
 * \param B Array of matrices 
 * \param f Bilinear map
 */

static inline void _mzd_ptr_apply_blm(mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f) {
	if (f->f!=NULL)
		_mzd_ptr_apply_blm_djb(X, A, B, f);
	else
		_mzd_ptr_apply_blm_mzd(X, A, B, f);
}


#endif //M4RIE_BLM_H
