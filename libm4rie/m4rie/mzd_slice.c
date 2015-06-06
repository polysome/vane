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

#include "mzd_slice.h"
#include "gf2x.h"
#include "mzd_ptr.h"
#include "m4ri_functions.h"

mzd_slice_t *mzd_slice_mul_scalar(mzd_slice_t *C, const word a, const mzd_slice_t *B) {
	if(C == NULL)
		C = mzd_slice_init(B->finite_field, B->nrows, B->ncols);
	else
		mzd_slice_set_ui(C, 0);
	assert( (C->finite_field == B->finite_field) && (((C->nrows ^ B->nrows) | (C->ncols ^ B->ncols)) == 0));

	const gf2e *ff = B->finite_field;

	for(int i=0; i<ff->degree; i++) {
		if(a&(1<<i)) {
			for(int j=0; j<B->depth; j++)
				_mzd_ptr_add_modred(ff, B->x[j], C->x, i+j);
		}
	}
	return C;
}

mzd_slice_t *mzd_slice_addmul_scalar(mzd_slice_t *C, const word a, const mzd_slice_t *B) {
	assert( (C->finite_field == B->finite_field) && (((C->nrows ^ B->nrows) | (C->ncols ^ B->ncols)) == 0));

	const gf2e *ff = B->finite_field;

	for(int i=0; i<ff->degree; i++) {
		if(a&(1<<i)) {
			for(int j=0; j<B->depth; j++)
				_mzd_ptr_add_modred(ff, B->x[j], C->x, i+j);
		}
	}
	return C;
}

void mzd_slice_set_ui(mzd_slice_t *A, word value) {
	for(int i=0; i<A->depth; i++) {
		mzd_set_ui(A->x[i], (value>>i)&1);
	}
}

void mzd_slice_print(const mzd_slice_t *A) {
	char formatstr[10];
	int width = gf2e_degree_to_w(A->finite_field)/4;
	if (gf2e_degree_to_w(A->finite_field)%4)
		width += 1;
	sprintf(formatstr,"%%%dx",width);

	for (rci_t i=0; i < A->nrows; ++i) {
		printf("[");
		for (rci_t j=0; j < A->ncols; j++) {
			word tmp = mzd_slice_read_elem(A,i,j);
			printf(formatstr,(int)tmp);
			if(j<A->ncols-1)
				printf(" ");
		}
		printf("]\n");
	}
}

mzd_slice_t *_mzd_slice_addmul_naive(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
	if (C == NULL)
		C = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

	const unsigned int e = A->finite_field->degree;

	mzd_t *t0 = mzd_init(A->nrows, B->ncols);

	for(unsigned int i=0; i<e; i++) {
		for(unsigned int j=0; j<e; j++) {
			mzd_mul(t0, A->x[i], B->x[j], 0);
			_mzd_ptr_add_modred(A->finite_field, t0, C->x, i+j);
		}
	}
	mzd_free(t0);
	return C;
}

