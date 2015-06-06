/**
 * \brief inline template for TRSM routines 
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 *
 * \note We want to keep this library in C, hence we cannot use of C++
 * templates.
 */

void _matrix_trsm_lower_left(const matrix_t *L, matrix_t *B, const rci_t cutoff) {
	assert((L->finite_field == B->finite_field) && (L->nrows == L->ncols) && (B->nrows == L->ncols));

	if (L->nrows <= cutoff || B->ncols <= cutoff) {
		matrix_trsm_lower_left_newton_john(L,B);
		return;
	}

	/**
	  \verbatim  
	  |\           ______
	  | \         |      |
	  |  \        |  B0  |
	  |L00\       |      |
	  |____\      |______|
	  |    |\     |      |
	  |    | \    |      |
	  |    |  \   |  B1  |
	  |L10 |L11\  |      |
	  |____|____\ |______|
	  \endverbatim 
	 */

	rci_t c = L->nrows/2;
	c = MAX((c - c%m4ri_radix),m4ri_radix);

	matrix_t *B0  = matrix_init_window(B, 0, 0, c, B->ncols);
	matrix_t *B1  = matrix_init_window(B, c, 0, B->nrows, B->ncols);
	const matrix_t *L00 = (const matrix_t*)matrix_init_window((matrix_t*)L, 0, 0, c, c);
	const matrix_t *L10 = (const matrix_t*)matrix_init_window((matrix_t*)L, c, 0, B->nrows, c);
	const matrix_t *L11 = (const matrix_t*)matrix_init_window((matrix_t*)L, c, c, B->nrows, B->nrows);

	_matrix_trsm_lower_left(L00, B0, cutoff);
	matrix_addmul(B1, L10, B0);
	_matrix_trsm_lower_left(L11, B1, cutoff);

	matrix_free_window(B0);
	matrix_free_window(B1);
	matrix_free_window((matrix_t*)L00);
	matrix_free_window((matrix_t*)L10);
	matrix_free_window((matrix_t*)L11);
}


void _matrix_trsm_upper_left(matrix_t const *U, matrix_t *B, const rci_t cutoff) {
	assert((U->finite_field == B->finite_field) && (U->nrows == U->ncols) && (B->nrows == U->ncols));

	if (U->nrows <= cutoff || B->ncols <= cutoff) {
		matrix_trsm_upper_left_newton_john(U,B);
		return;
	}
	/**
	  \verbatim
	  __________   ______
	  \ U00|    | |      |
	  \   |U01 | |      |
	  \  |    | |  B0  |
	  \ |    | |      |
	  \|____| |______|
	  \    | |      |
	  \U11| |      |
	  \  | |  B1  |
	  \ | |      |
	  \| |______|
	  \endverbatim 
	 */

	rci_t c = U->nrows/2;
	c = MAX((c - c%m4ri_radix),m4ri_radix);

	matrix_t *B0  = matrix_init_window(B,  0,  0, c, B->ncols);
	matrix_t *B1  = matrix_init_window(B,  c, 0, B->nrows, B->ncols);
	const matrix_t *U00 = (const matrix_t *)matrix_init_window(U,  0,  0, c, c);
	const matrix_t *U01 = (const matrix_t *)matrix_init_window(U,  0, c, c, B->nrows);
	const matrix_t *U11 = (const matrix_t *)matrix_init_window(U, c, c, B->nrows, B->nrows);

	_matrix_trsm_upper_left(U11, B1, cutoff);
	matrix_addmul(B0, U01, B1);
	_matrix_trsm_upper_left(U00, B0, cutoff);

	matrix_free_window(B0);
	matrix_free_window(B1);    
	matrix_free_window((matrix_t*)U00);
	matrix_free_window((matrix_t*)U01);
	matrix_free_window((matrix_t*)U11);
}


