#include "trsm.h"
#include "newton_john.h"
#include "conversion.h"

void mzed_trsm_upper_left_naive(const mzed_t *U, mzed_t *B) {
	assert(U->finite_field == B->finite_field);
	assert(U->nrows == U->ncols);
	assert(B->nrows == U->ncols);

	const gf2e *ff = U->finite_field;
	for(int i=B->nrows-1; i>=0; i--) {
		for(rci_t k=i+1; k<B->nrows; k++) {
			mzed_add_multiple_of_row(B, i, B, k, mzed_read_elem(U, i, k), 0);
		}
		mzed_rescale_row(B, i, 0, gf2e_inv(ff, mzed_read_elem(U, i, i)));
	}
}

void mzed_trsm_lower_left_naive(const mzed_t *L, mzed_t *B) {
	assert(L->finite_field == B->finite_field);
	assert(L->nrows == L->ncols);
	assert(B->nrows == L->ncols);

	const gf2e *ff = L->finite_field;
	for(rci_t i=0; i<B->nrows; i++) {
		for(rci_t k=0; k<i; k++) {
			mzed_add_multiple_of_row(B, i, B, k, mzed_read_elem(L, i, k), 0);
		}
		mzed_rescale_row(B, i, 0, gf2e_inv(ff, mzed_read_elem(L, i, i)));
	}
}

void mzd_slice_trsm_upper_left_naive(const mzd_slice_t *U, mzd_slice_t *B) {
	assert(U->finite_field == B->finite_field);
	assert(U->nrows == U->ncols);
	assert(B->nrows == U->ncols);

	const mzed_t *Ue = mzed_cling(NULL, U);
	mzed_t *Be       = mzed_cling(NULL, B);
	mzed_trsm_upper_left_naive(Ue, Be);

	mzed_slice(B, Be);
	mzed_free((mzed_t*)Ue);
	mzed_free(Be);
}

void mzd_slice_trsm_lower_left_naive(const mzd_slice_t *L, mzd_slice_t *B) {
	assert(L->finite_field == B->finite_field);
	assert(L->nrows == L->ncols);
	assert(B->nrows == L->ncols);

	const mzed_t *Le = mzed_cling(NULL, L);
	mzed_t *Be = mzed_cling(NULL, B);

	mzed_trsm_lower_left_naive(Le, Be);

	mzed_slice(B, Be);
	mzed_free((mzed_t*)Le);
	mzed_free(Be);
}

#include "mzed_intro.inl"
#include "trsm.inl"
#include "mzed_outro.inl"

#include "mzd_slice_intro.inl"
#include "trsm.inl"
#include "mzd_slice_outro.inl"

