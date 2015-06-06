/**
 * \file newton_john.h
 *
 * \brief Newton-John table based algorithms
 *
 * \note These tables were formally known as Travolta tables.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_NEWTON_JOHN_H
#define M4RIE_NEWTON_JOHN_H

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

#include <m4rie/gf2e.h>
#include <m4rie/mzed.h>
#include <m4rie/mzd_slice.h>

/**
 * \brief Newton-John table
 */

typedef struct {
	rci_t *L;  /**< A map such that L[a] points to the row where the first entry is a. */
	mzed_t *M; /**< Table of length \e with multiples of the input s.t. \f$a^i\f$ is the first entry of row \f$i\f$. */
	mzed_t *T; /**< Actual table of length \f$2^e\f$ of all linear combinations of T. */
} njt_mzed_t;

/**
 * \brief Allocate Newton-John table of dimension gf2e::degree<<1 * ncols.
 *
 * \param ff Finite field.
 * \param ncols Integer > 0.
 */

njt_mzed_t *njt_mzed_init(const gf2e *ff, const rci_t ncols);

/**
 * \brief Free Newton-John table
 *
 * \param t Table
 */

void njt_mzed_free(njt_mzed_t *t);

/**
 * \brief Construct Newton-John table T for row r of A, and element A[r,c].
 *
 * \param T Preallocated Newton-John table or NULL.
 * \param A Matrix.
 * \param r Row index.
 * \param c Column index.
 */

njt_mzed_t * mzed_make_table(njt_mzed_t *T, const mzed_t *A, const rci_t r, const rci_t c);

/**
 * \brief \f$C = A \cdot B\f$ using Newton-John tables.
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa mzed_mul _mzed_mul_newton_john0()
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_mul_newton_john(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$C = C + A \cdot B\f$ using Newton-John tables.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzed_mul_newton_john() mzed_mul()
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_addmul_newton_john(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$C = C + A \cdot B\f$ using Newton-John tables.
 *
 * This is a simple implementation for clarity of presentation. Do not
 * call, it is slow.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa mzed_mul_newton_john() mzed_mul()
 *
 * \ingroup Multiplication
 */

mzed_t *_mzed_mul_newton_john0(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$C = C + A \cdot B\f$ using Newton-John tables.
 *
 * This is an optimised implementation.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa mzed_mul()
 *
 * \ingroup Multiplication
 */

mzed_t *_mzed_mul_newton_john(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Reduce matrix A to row echelon form using Gauss-Newton-John
 * elimination.
 *
 * \param A Matrix to be reduced.
 * \param full If set to true, the reduced row echelon form will be
 * computed.
 *
 * \ingroup Echelon
 */

rci_t mzed_echelonize_newton_john(mzed_t *A, int full);

/**
 * \brief Invert the matrix A using Gauss-Newton-John elimination.
 *
 * \param B Preallocated space for inversion matrix, may be NULL for
 * automatic creation.
 * \param A Matrix to be inverted.
 */

mzed_t *mzed_invert_newton_john(mzed_t *B, const mzed_t *A);

/**
 * \brief \f$B = L^{-1} \cdot B\f$ using Newton-John tables.
 *
 * \param L Lower-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

void mzed_trsm_lower_left_newton_john(const mzed_t *L, mzed_t *B);

/**
 * \brief \f$B = L^{-1} \cdot B\f$ using Newton-John tables.
 *
 * \param L Lower-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

void mzd_slice_trsm_lower_left_newton_john(const mzd_slice_t *L, mzd_slice_t *B);

/**
 * \brief \f$B = U^{-1} \cdot B\f$ using Newton-John tables.
 *
 * \param U Upper-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

void mzed_trsm_upper_left_newton_john(const mzed_t *U, mzed_t *B);

/**
 * \brief \f$B = U^{-1} \cdot B\f$ using Newton-John tables.
 *
 * \param U Upper-triangular matrix (other entries are ignored).
 * \param B Matrix.
 *
 * \ingroup Triangular
 */

void mzd_slice_trsm_upper_left_newton_john(const mzd_slice_t *U, mzd_slice_t *B);

/**
 * \brief PLE decomposition: \f$L \cdot E = P\cdot A\f$ using Newton-John tables.
 *
 * \ingroup PLE
 */

rci_t mzed_ple_newton_john(mzed_t *A, mzp_t *P, mzp_t *Q);

/**
 * \brief The function looks up 6 entries from position i,startcol in
 * each row and adds the appropriate row from T to the row i.
 *
 * This process is iterated for i from startrow to stoprow
 * (exclusive).
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T Newton-John table
 *
 * \ingroup RowOperations
 */

static inline void mzed_process_rows(mzed_t *M, const rci_t startrow, const rci_t endrow, rci_t startcol, const njt_mzed_t *T) {
	mzd_process_rows(M->x, startrow, endrow, startcol*M->w, M->w, T->T->x, T->L);
}

/**
 * \brief Same as mzed_process_rows but works with two Newton-John tables
 * in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Newton-John table
 * \param T1 Newton-John table
 *
 * \ingroup RowOperations
 */

static inline void mzed_process_rows2(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol, 
		const njt_mzed_t *T0, const njt_mzed_t *T1) {
	mzd_process_rows2(M->x, startrow, endrow, startcol*M->w, 2*M->w, T0->T->x, T0->L, T1->T->x, T1->L);
}

/**
 * \brief Same as mzed_process_rows but works with three Newton-John
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Newton-John table
 * \param T1 Newton-John table
 * \param T2 Newton-John table
 *
 * \ingroup RowOperations
 */

static inline void mzed_process_rows3(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol,
		const njt_mzed_t *T0, const njt_mzed_t *T1, const njt_mzed_t *T2) {
	mzd_process_rows3(M->x, startrow, endrow, startcol*M->w, 3*M->w, T0->T->x, T0->L, T1->T->x, T1->L, T2->T->x, T2->L);
}

/**
 * \brief Same as mzed_process_rows but works with four Newton-John
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Newton-John table
 * \param T1 Newton-John table
 * \param T2 Newton-John table
 * \param T3 Newton-John table
 *
 * \ingroup RowOperations
 */

static inline void mzed_process_rows4(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol,
		const njt_mzed_t *T0, const njt_mzed_t *T1, const njt_mzed_t *T2, const njt_mzed_t *T3) {
	mzd_process_rows4(M->x, startrow, endrow, startcol*M->w, 4*M->w, T0->T->x, T0->L, T1->T->x, T1->L, T2->T->x, T2->L, T3->T->x, T3->L);
}

/**
 * \brief Same as mzed_process_rows but works with five Newton-John
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Newton-John table
 * \param T1 Newton-John table
 * \param T2 Newton-John table
 * \param T3 Newton-John table
 * \param T4 Newton-John table
 *
 * \ingroup RowOperations
 */

static inline void mzed_process_rows5(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol,
		const njt_mzed_t *T0, const njt_mzed_t *T1, const njt_mzed_t *T2, const njt_mzed_t *T3, const njt_mzed_t *T4) {
	mzd_process_rows5(M->x, startrow, endrow, startcol*M->w, 5*M->w, T0->T->x, T0->L, T1->T->x, T1->L, T2->T->x, T2->L, T3->T->x, T3->L, T4->T->x, T4->L);
}


/**
 * \brief Same as mzed_process_rows but works with six Newton-John tables
 * in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Newton-John table
 * \param T1 Newton-John table
 * \param T2 Newton-John table
 * \param T3 Newton-John table
 * \param T4 Newton-John table
 * \param T5 Newton-John table
 *
 * \ingroup RowOperations
 */

static inline void mzed_process_rows6(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol,
		const njt_mzed_t *T0, const njt_mzed_t *T1, const njt_mzed_t *T2,
		const njt_mzed_t *T3, const njt_mzed_t *T4, const njt_mzed_t *T5) {
	mzd_process_rows6(M->x, startrow, endrow, startcol*M->w, 6*M->w, T0->T->x, T0->L, T1->T->x, T1->L, T2->T->x, T2->L, T3->T->x, T3->L, T4->T->x, T4->L, T5->T->x, T5->L);
}


#endif //M4RIE_NEWTON_JOHN_H
