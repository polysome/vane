/**
 * \file strassen.h
 * \brief Strassen-Winograd multiplication for mzed_t
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_STRASSEN_H
#define M4RIE_STRASSEN_H

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

/**
 * \brief \f$ C = A \cdot B \f$ using Strassen-Winograd.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Newton-John table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix, may be NULL for allocation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 * \param cutoff Crossover to basecase dimension > 64 or 0 for heuristic choice
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_mul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff);

/**
 * \brief \f$ C = C + A \cdot B \f$ using Strassen-Winograd.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Newton-John table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix, may be NULL for allocation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 * \param cutoff Crossover to basecase dimension > 64 or 0 for heuristic choice.
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_addmul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff);

/**
 * \brief \f$ C = A \cdot B \f$ using Strassen-Winograd.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Newton-John table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 * \param cutoff Crossover to basecase dimension > 64
 *
 * \ingroup Multiplication
 *
 */

mzed_t *_mzed_mul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff);

/**
 * \brief \f$ C = A \cdot B \f$ using Strassen-Winograd.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Newton-John table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 * \param cutoff Crossover to basecase dimension > 64
 *
 * \ingroup Multiplication
 */

mzed_t *_mzed_addmul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff);

/**
 * \brief Return heurstic choice for crossover parameter for Strassen-Winograd multiplication given A, B and C.
 *
 * \param C Matrix (ignored)
 * \param A Matrix
 * \param B Martix (ignored)
 *
 * \ingroup Multiplication
 */

rci_t _mzed_strassen_cutoff(const mzed_t *C, const mzed_t *A, const mzed_t *B);

#endif //M4RIE_STRASSEN_H
