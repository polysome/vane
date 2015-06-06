#include "mzd_ptr.h"
#include "m4ri_functions.h"


void djb_apply_mzd_ptr(djb_t *m, mzd_t **W, const mzd_t **V) {
	int *iszero = m4ri_mm_malloc(sizeof(int)*m->nrows);
	for(int i = 0; i < m->nrows; ++i) 
		iszero[i] = 1;

	rci_t i = m->length;
	while (i > 0) {
		--i;
		if (iszero[m->target[i]]) {
			if (m->srctyp[i] == source_source) {
				mzd_copy(W[m->target[i]], V[m->source[i]]);
			} else {
				mzd_copy(W[m->target[i]], W[m->source[i]]);
			}
			iszero[m->target[i]] = 0;
		} else {
			if (m->srctyp[i] == source_source) {
				mzd_add(W[m->target[i]], W[m->target[i]], V[m->source[i]]);
			} else {
				mzd_add(W[m->target[i]], W[m->target[i]], W[m->source[i]]);
			}
		}
	}
	m4ri_mm_free(iszero);
}
