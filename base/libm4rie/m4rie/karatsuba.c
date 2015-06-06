#include "mzd_ptr.h"
#include "m4ri_functions.h"

void _mzd_ptr_addmul_karatsuba2(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];

	mzd_t *t0 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t1 = mzd_init(b0->nrows, b0->ncols);

	mzd_addmul(X[1], mzd_add(t0, a0, a1), mzd_add(t1, b0, b1), 0); /* + (a0+a1)(b0+b1)X */

	mzd_free(t0);
	mzd_free(t1);

	t0 = mzd_init(a0->nrows, b0->ncols);
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 2,  0, 1); /* + a0b0(1-X) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 2,  1, 2); /* + a1b1(X+X^2) */
	mzd_free(t0);
}

void _mzd_ptr_addmul_karatsuba3(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_add(t1, a0, a1), mzd_add(t2, b0, b1), 0), X, 1,   1); /* + (a0+a1)(b0+b1)X */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_add(t1, a0, a2), mzd_add(t2, b0, b2), 0), X, 1,   2); /* + (a0+a2)(b0+b2)X^2 */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_add(t1, a1, a2), mzd_add(t2, b1, b2), 0), X, 1,   3); /* + (a1+a2)(b1+b2)X^3 */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 3,   0, 1, 2); /* + a0b0(1-X-X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 3,   1, 2, 3); /* + a1b1(X+X^2-X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a2, b2, 0), X, 3,   2, 3, 4); /* + a2b2(-X^2-X^3+X^4) */

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba4(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	/* (a0 + a1 + a2 + a3)*(b0 + b1 + b2 + b3)*X^3 */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a2, a3),
				mzd_sum(t2, 4, b0, b1, b2, b3), 0), X, 1,   3);
	/* (a0 + a1)*(b0 + b1)*(X   + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 2,   1, 3);
	/* (a0 + a2)*(b0 + b2)*(X^2 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a2),
				mzd_sum(t2, 2, b0, b2), 0), X, 2,   2, 3);
	/* (a1 + a3)*(b1 + b3)*(X^3 + X^4) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a1, a3),
				mzd_sum(t2, 2, b1, b3), 0), X, 2,   3, 4);
	/* (a2 + a3)*(b2 + b3)*(X^3 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a2, a3),
				mzd_sum(t2, 2, b2, b3), 0), X, 2,   3, 5);
	/* (a0*b0)*(1 + X   + X^2 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 4,    0, 1, 2, 3);
	/* (a1*b1)*(X + X^2 + X^3 + X^4) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 4,    1, 2, 3, 4);
	/* (a2*b2)*(X^2 + X^3 + X^4 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a2, b2, 0), X, 4,    2, 3, 4, 5);
	/* (a3*b3)*(X^3 + X^4 + X^5 + X^6) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a3, b3, 0), X, 4,    3, 4, 5, 6);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba5(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];
	const mzd_t *a4 = A[4];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];
	const mzd_t *b4 = B[4];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);

	/* a0b0(X^6 + X^5 + X^4 + X^2 + X + 1)*/
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 6,   6, 5, 4, 2, 1, 0);
	/* a1b1(X^4 + X)*/
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 2,   4, 1);
	/* a3b3(X^7 + X^4) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a3, b3, 0), X, 2,   7, 4);
	/* (a4b4)(X^8 + X^7 + X^6 + X^4 + X^3 + X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a4, b4, 0), X, 6,   8, 7, 6, 4, 3, 2);

	mzd_t *t1 = mzd_init(a0->nrows, a1->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b1->ncols);

	/* (a0+a4)(b0+b4)(X^6 + X^5 + X^3 + X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a4),
				mzd_sum(t2, 2, b0, b4), 0), X, 4,   6, 5, 3, 2);

	/* (a0+a1)(b0+b1)(X^5 + X^4 + X^2 + X) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 4,   5, 4, 2, 1);

	/* (a3+a4)(b3+b4)(X^7 + X^6 + X^4 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a3, a4),
				mzd_sum(t2, 2, b3, b4), 0), X, 4,   7, 6, 4, 3);

	/* (a1+a2+a4)(a1+a2+a4)(X^4 + X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a1, a2, a4),
				mzd_sum(t2, 3, b1, b2, b4), 0), X, 2,   4, 2);

	/* (a0+a2+a3)(b0+b2+b3)(X^6 + X^4) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a0, a2, a3),
				mzd_sum(t2, 3, b0, b2, b3), 0), X, 2,   6, 4);

	/* (a0+a1+a3+a4)(b0+b1+b3+b4)(X^5 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a3, a4),
				mzd_sum(t2, 4, b0, b1, b3, b4), 0), X, 2,   5, 3);

	/* (a0+a1+a2+a4)(b0+b1+b2+b4)(X^5 + X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a2, a4),
				mzd_sum(t2, 4, b0, b1, b2, b4), 0), X, 2,   5, 2);

	/* (a0+a2+a3+a4)(b0+b2+b3+b4)(X^6 + X^3)*/
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a2, a3, a4),
				mzd_sum(t2, 4, b0, b2, b3, b4), 0), X, 2,   6, 3);

	/* (a0+a1+a2+a3+a4)(b0+b1+b2+b3+b4)(X^5 + X^4 + X^3)*/
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a2, a3, a4),
				mzd_sum(t2, 5, b0, b1, b2, b3, b4), 0), X, 3,   5, 4, 3);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba6(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];
	const mzd_t *a4 = A[4];
	const mzd_t *a5 = A[5];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];
	const mzd_t *b4 = B[4];
	const mzd_t *b5 = B[5];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a1->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b1->ncols);

	/* a0b0 (X^6 + X^5 + X + 1) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 4,   6, 5, 1, 0);
	/* a1b1 (X^7 + X^6 + X^5 + X^4 + X^3 + X) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 6,   7, 6, 5, 4, 3, 1);
	/* a4b4 (X^9 + X^7 + X^5 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a4, b4, 0), X, 4,   9, 7, 5, 3);
	/* a5b5 (X^10 + X^9 + X^6 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a5, b5, 0), X, 4,   10, 9, 6, 5);

	/* (a4 + a5)(b4 + b5) (X^9 + X^8 + X^4+ X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a4, a5),
				mzd_sum(t2, 2, b4, b5), 0), X, 4,   9, 8, 4, 3);
	/* (a0 + a1)(b0 + b1) (X^7 + X^4 + X^2 + X) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 4,   7, 4, 2, 1);
	/* (a3 + a4)(b3 + b4)(X^8 + X^7 + X^6 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a3, a4),
				mzd_sum(t2, 2, b3, b4), 0), X, 4,   8, 7, 6, 3);
	/* (a1 + a2)(b1 + b2) (X^7 + X^6 + X^3 + X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a1, a2),
				mzd_sum(t2, 2, b1, b2), 0), X, 4,   7, 6, 3, 2);
	/* (a1 + a4)(b1 + b4) (X^4 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a1, a4),
				mzd_sum(t2, 2, b1, b4), 0), X, 2,   4, 5);
	/* (a2 + a3)(b2 + b3) (X^7 + X^6 + X^4 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a2, a3),
				mzd_sum(t2, 2, b2, b3), 0), X, 4,   7, 6, 4, 3);
	/* (a3 + a4 + a5)(b3 + b4 + b5) (X^8 + X^6 + X^4 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a3, a4, a5),
				mzd_sum(t2, 3, b3, b4, b5), 0), X, 4,   8, 6, 4, 3);
	/* (a0 + a1 + a2)(b0 + b1 + b2) (X^7 + X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a0, a1, a2),
				mzd_sum(t2, 3, b0, b1, b2), 0), X, 2,   7, 2);
	/* (a0 + a3 + a5)(b0 + b3 + b5) (X^7 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a0, a3, a5),
				mzd_sum(t2, 3, b0, b3, b5), 0), X, 2,   7, 5);
	/* (a0 + a2 + a5)(b0 + b2 + b5) (X^6 + X^5 + X^4 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a0, a2, a5),
				mzd_sum(t2, 3, b0, b2, b5), 0), X, 4,   6, 5, 4, 3);
	/* (a0 + a2 + a3 + a5)(b0 + b2 + b3 + b5) (X^7 + X^5 + X^4 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a2, a3, a5),
				mzd_sum(t2, 4, b0, b2, b3, b5), 0), X, 4,   7, 5, 4, 3);
	/* (a0 + a1 + a3 + a4)(b0 + b1 + b3 + b4) (X^6 + X^4) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a3, a4),
				mzd_sum(t2, 4, b0, b1, b3, b4), 0), X, 2,   6, 4);
	/* (a0 + a1 + a2 + a3 + a4 + a5)(b0 + b1 + b2 + b3 + b4 + b5) X^6 */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a2, a3, a4, a5),
				mzd_sum(t2, 6, b0, b1, b2, b3, b4, b5), 0), X, 1,   6);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba7(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 22 multiplications and three temporaries */
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];
	const mzd_t *a4 = A[4];
	const mzd_t *a5 = A[5];
	const mzd_t *a6 = A[6];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];
	const mzd_t *b4 = B[4];
	const mzd_t *b5 = B[5];
	const mzd_t *b6 = B[6];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a1->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b1->ncols);

	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6)(b0 + b1 + b2 + b3 + b4 + b5 + b6)(X^7 + X^6 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a2, a3, a4, a5, a6),
				mzd_sum(t2, 7, b0, b1, b2, b3, b4, b5, b6), 0),
			X, 3,   7, 6, 5);

	/* (a1 + a2 + a3 + a5 + a6)(b1 + b2 + b3 + b5 + b6)(X^9 + X^6) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a1, a2, a3, a5, a6),
				mzd_sum(t2, 5, b1, b2, b3, b5, b6), 0),
			X, 2,   9, 6);

	/* (a0 + a1 + a3 + a4 + a5)(b0 + b1 + b3 + b4 + b5)(X^6 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a3, a4, a5),
				mzd_sum(t2, 5, b0, b1, b3, b4, b5), 0),
			X, 2,   6, 3);

	/* (a0 + a2 + a3 + a4 + a6)(b0 + b2 + b3 + b4 + b6)(X^9 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a3, a4, a6),
				mzd_sum(t2, 5, b0, b2, b3, b4, b6), 0),
			X, 2,   9, 3);

	/* (a0 + a2 + a3 + a5 + a6)(b0 + b2 + b3 + b5 + b6)(X^7 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a3, a5, a6),
				mzd_sum(t2, 5, b0, b2, b3, b5, b6), 0),
			X, 2,   7, 3);

	/* (a0 + a1 + a3 + a4 + a6)(b0 + b1 + b3 + b4 + b6)(X^9 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a3, a4, a6),
				mzd_sum(t2, 5, b0, b1, b3, b4, b6), 0),
			X, 2,   9, 5);

	/* (a1 + a2 + a4 + a5)(b1 + b2 + b4 + b5)(X^9 + X^7 + X^5 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a1, a2, a4, a5),
				mzd_sum(t2, 4, b1, b2, b4, b5), 0),
			X, 4,   9, 7, 5, 3);

	/* (a0 + a1)(b0 + b1)(X^9 + X^7 + X^3 + X) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0),
			X, 4,   9, 7, 3, 1);

	/* (a0 + a2)(b0 + b2)(X^9 + X^6 + X^5 + X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a2),
				mzd_sum(t2, 2, b0, b2), 0),
			X, 4,   9, 6, 5, 2);

	/* (a0 + a4)(b0 + b4)(X^7 + X^4) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a4),
				mzd_sum(t2, 2, b0, b4), 0),
			X, 2,   7, 4);

	/* (a1 + a3)(b1 + b3)(X^7 + X^6 + X^4 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a1, a3),
				mzd_sum(t2, 2, b1, b3), 0),
			X, 4,   7, 6, 4, 3);

	/* (a2 + a6)(b2 + b6)(X^8 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a2, a6),
				mzd_sum(t2, 2, b2, b6), 0),
			X, 2,   8, 5);

	/* (a3 + a5)(b3 + b5)(X^9 + X^8 + X^6 + X^5) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a3, a5),
				mzd_sum(t2, 2, b3, b5), 0),
			X, 4,   9, 8, 6, 5);

	/* (a4 + a6)(b4 + b6)(X^10 + X^7 + X^6 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a4, a6),
				mzd_sum(t2, 2, b4, b6), 0),
			X, 4,  10, 7, 6, 3);

	/* (a5 + a6)(b5 + b6)(X^11 + X^9 + X^5 + X^3) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a5, a6),
				mzd_sum(t2, 2, b5, b6), 0),
			X, 4,  11, 9, 5, 3);

	/* a0b0(X^6 + X^5 + X^4 + X^2 + X + 1) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 6,   6, 5, 4, 2, 1, 0);
	/* a1b1(X^5 + X^4 + X^2 + X) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 4,   5, 4, 2, 1);
	/* a2b2(X^8 + X^7 + X^6 + X^4 + X^3 + X^2) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a2, b2, 0), X, 6,   8, 7, 6, 4, 3, 2);
	/* a3b3(X^8 + X^7 + X^5 + X^4) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a3, b3, 0), X, 4,   8, 7, 5, 4);
	/* a4b4(X^10 + X^9 + X^8 + X^6 + X^5 + X^4)*/
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a4, b4, 0), X, 6,  10, 9, 8, 6, 5, 4);;
	/* a5b5(X^11 + X^10 + X^8 + X^7) */
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a5, b5, 0), X, 4,  11, 10, 8, 7);
	/* a6b6(X^12 + X^11 + X^10 + X^8 + X^7 + X^6)*/
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a6, b6, 0), X, 6,  12, 11, 10, 8, 7, 6);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);

}

void _mzd_ptr_addmul_karatsuba8(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 27 multiplications and three temporaries (not optimal!) */
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];
	const mzd_t *a4 = A[4];
	const mzd_t *a5 = A[5];
	const mzd_t *a6 = A[6];
	const mzd_t *a7 = A[7];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];
	const mzd_t *b4 = B[4];
	const mzd_t *b5 = B[5];
	const mzd_t *b6 = B[6];
	const mzd_t *b7 = B[7];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	// a0*b0 * (X^0 + X^1 + X^2 + X^3 + X^4 + X^5 + X^6 + X^7) 
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 8,  0, 1, 2, 3, 4, 5, 6, 7);
	// a1*b1 * (X^1 + X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^8) 
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 8,  1, 2, 3, 4, 5, 6, 7, 8);
	// a2*b2 * (X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X^9) 
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a2, b2, 0), X, 8,  2, 3, 4, 5, 6, 7, 8, 9);
	// a3*b3 * (X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X^10) 
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a3, b3, 0), X, 8,  3, 4, 5, 6, 7, 8, 9, 10);
	// a4*b4 * (X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X^10 + X^11) 
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a4, b4, 0), X, 8,  4, 5, 6, 7, 8, 9, 10, 11);
	// a5*b5 * (X^5 + X^6 + X^7 + X^8 + X^9 + X^10 + X^11 + X^12) 
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a5, b5, 0), X, 8,  5, 6, 7, 8, 9, 10, 11, 12);
	// a6*b6 * (X^6 + X^7 + X^8 + X^9 + X^10 + X^11 + X^12 + X^13) 
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a6, b6, 0), X, 8,  6, 7, 8, 9, 10, 11, 12, 13);
	// a7*b7 * (X^7 + X^8 + X^9 + X^10 + X^11 + X^12 + X^13 + X^14) 
	_mzd_ptr_add_to_all(ff, mzd_mul(t0, a7, b7, 0), X, 8,  7, 8, 9, 10, 11, 12, 13, 14);

	// (b0 + b2)*(a0 + a2) * (X^2 + X^3 + X^6 + X^7)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a2),
				mzd_sum(t2, 2, b0, b2), 0), X, 4,  2, 3, 6, 7);
	// (b2 + b3)*(a2 + a3) * (X^3 + X^5 + X^7 + X^9)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a2, a3),
				mzd_sum(t2, 2, b2, b3), 0), X, 4,  3, 5, 7, 9);
	// (b2 + b6)*(a2 + a6) * (X^6 + X^7 + X^8 + X^9)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a2, a6),
				mzd_sum(t2, 2, b2, b6), 0), X, 4,  6, 7, 8, 9);
	// (b4 + b5)*(a4 + a5) * (X^5 + X^7 + X^9 + X^11)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a4, a5),
				mzd_sum(t2, 2, b4, b5), 0), X, 4,  5, 7, 9, 11);
	// (b6 + b7)*(a6 + a7) * (X^7 + X^9 + X^11 + X^13)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a6, a7),
				mzd_sum(t2, 2, b6, b7), 0), X, 4,  7, 9, 11, 13);
	// (b0 + b1 + b2 + b3)*(a0 + a1 + a2 + a3) * (X^3 + X^7)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a2, a3),
				mzd_sum(t2, 4, b0, b1, b2, b3), 0), X, 2,  3, 7);
	// (b5 + b7)*(a5 + a7) * (X^7 + X^8 + X^11 + X^12)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a5, a7),
				mzd_sum(t2, 2, b5, b7), 0), X, 4,  7, 8, 11, 12);
	// (b0 + b4)*(a0 + a4) * (X^4 + X^5 + X^6 + X^7)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a4),
				mzd_sum(t2, 2, b0, b4), 0), X, 4,  4, 5, 6, 7);
	// (b1 + b3)*(a1 + a3) * (X^3 + X^4 + X^7 + X^8)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a1, a3),
				mzd_sum(t2, 2, b1, b3), 0), X, 4,  3, 4, 7, 8);
	// (b1 + b3 + b5 + b7)*(a1 + a3 + a5 + a7) * (X^7 + X^8)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a1, a3, a5, a7),
				mzd_sum(t2, 4, b1, b3, b5, b7), 0), X, 2,  7, 8);
	// (b0 + b1)*(a0 + a1) * (X^1 + X^3 + X^5 + X^7)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 4,  1, 3, 5, 7);
	// (b4 + b6)*(a4 + a6) * (X^6 + X^7 + X^10 + X^11)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a4, a6),
				mzd_sum(t2, 2, b4, b6), 0), X, 4,  6, 7, 10, 11);
	// (b1 + b5)*(a1 + a5) * (X^5 + X^6 + X^7 + X^8)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a1, a5),
				mzd_sum(t2, 2, b1, b5), 0), X, 4,  5, 6, 7, 8);
	// (b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7)*(a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7) * (X^7)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a2, a3, a4, a5, a6, a7),
				mzd_sum(t2, 8, b0, b1, b2, b3, b4, b5, b6, b7), 0), X, 1,  7);
	// (b0 + b2 + b4 + b6)*(a0 + a2 + a4 + a6) * (X^6 + X^7)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a2, a4, a6),
				mzd_sum(t2, 4, b0, b2, b4, b6), 0), X, 2,  6, 7);
	// (b0 + b1 + b4 + b5)*(a0 + a1 + a4 + a5) * (X^5 + X^7)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a4, a5),
				mzd_sum(t2, 4, b0, b1, b4, b5), 0), X, 2,  5, 7);
	// (b2 + b3 + b6 + b7)*(a2 + a3 + a6 + a7) * (X^7 + X^9)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a2, a3, a6, a7),
				mzd_sum(t2, 4, b2, b3, b6, b7), 0), X, 2,  7, 9);
	// (b3 + b7)*(a3 + a7) * (X^7 + X^8 + X^9 + X^10)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a3, a7),
				mzd_sum(t2, 2, b3, b7), 0), X, 4,  7, 8, 9, 10);
	// (b4 + b5 + b6 + b7)*(a4 + a5 + a6 + a7) * (X^7 + X^11)
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a4, a5, a6, a7),
				mzd_sum(t2, 4, b4, b5, b6, b7), 0), X, 2,  7, 11);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba9(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 31 multiplications and three temporaries (not optimal!) */

	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];
	const mzd_t *a4 = A[4];
	const mzd_t *a5 = A[5];
	const mzd_t *a6 = A[6];
	const mzd_t *a7 = A[7];
	const mzd_t *a8 = A[8];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];
	const mzd_t *b4 = B[4];
	const mzd_t *b5 = B[5];
	const mzd_t *b6 = B[6];
	const mzd_t *b7 = B[7];
	const mzd_t *b8 = B[8];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	/* (a8)(b8)*(X^2 + X^3 + X^5 + X0 + X2 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a8, b8, 0), X, 6,  2, 3, 5, 10, 12, 16);
	/* (a0)(b0)*(1 + X + X^3 + X^8 + X0 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a0, b0, 0), X, 6,  0, 1, 3, 8, 10, 14);
	/* (a1)(b1)*(X + X^2 + X^4 + X^9 + X1 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a1, b1, 0), X, 6,  1, 2, 4, 9, 11, 15);
	/* (a0 + a1)(b0 + b1)*(X + X^2 + X^4 + X^9 + X1 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 6,  1, 2, 4, 9, 11, 15);
	/* (a0 + a2 + a4 + a6 + a8)(b0 + b2 + b4 + b6 + b8)*(X^2 + X^5 + X^6 + X^7 + X^8 + X^9 + X2 + X3 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a4, a6, a8),
				mzd_sum(t2, 5, b0, b2, b4, b6, b8), 0), X, 10,  2, 5, 6, 7, 8, 9, 12, 13, 14, 15);
	/* (a1 + a3 + a5 + a7)(b1 + b3 + b5 + b7)*(X^2 + X^5 + X^6 + X^7 + X^8 + X^9 + X2 + X3 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a1, a3, a5, a7),
				mzd_sum(t2, 4, b1, b3, b5, b7), 0), X, 10,  2, 5, 6, 7, 8, 9, 12, 13, 14, 15);
	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8)(b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8)*(X^2 + X^3 + X^4 + X^6 + X^8 + X2 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a2, a3, a4, a5, a6, a7, a8),
				mzd_sum(t2, 9, b0, b1, b2, b3, b4, b5, b6, b7, b8), 0), X, 7,  2, 3, 4, 6, 8, 12, 14);
	/* (a0 + a2 + a3 + a5 + a6 + a8)(b0 + b2 + b3 + b5 + b6 + b8)*(X^2 + X^3 + X^4 + X^5 + X^6 + X^8 + X^9 + X0 + X1 + X2 + X3 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a3, a5, a6, a8),
				mzd_sum(t2, 6, b0, b2, b3, b5, b6, b8), 0), X, 12,  2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15);
	/* (a1 + a2 + a4 + a5 + a7 + a8)(b1 + b2 + b4 + b5 + b7 + b8)*(X^2 + X^4 + X^6 + X^7 + X^9 + X1 + X3 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a2, a4, a5, a7, a8),
				mzd_sum(t2, 6, b1, b2, b4, b5, b7, b8), 0), X, 8,  2, 4, 6, 7, 9, 11, 13, 14);
	/* (a0 + a1 + a3 + a4 + a6 + a7)(b0 + b1 + b3 + b4 + b6 + b7)*(X^3 + X^5 + X^7 + X^8 + X0 + X2 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a3, a4, a6, a7),
				mzd_sum(t2, 6, b0, b1, b3, b4, b6, b7), 0), X, 8,  3, 5, 7, 8, 10, 12, 14, 15);
	/* (a0 + a3 + a5 + a6 + a7)(b0 + b3 + b5 + b6 + b7)*(X^2 + X^3 + X0 + X2 + X3 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a3, a5, a6, a7),
				mzd_sum(t2, 5, b0, b3, b5, b6, b7), 0), X, 6,  2, 3, 10, 12, 13, 14);
	/* (a1 + a3 + a4 + a5 + a8)(b1 + b3 + b4 + b5 + b8)*(X^3 + X^4 + X1 + X3 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a3, a4, a5, a8),
				mzd_sum(t2, 5, b1, b3, b4, b5, b8), 0), X, 6,  3, 4, 11, 13, 14, 15);
	/* (a2 + a4 + a5 + a6)(b2 + b4 + b5 + b6)*(X^2 + X^3 + X^4 + X0 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a2, a4, a5, a6),
				mzd_sum(t2, 4, b2, b4, b5, b6), 0), X, 6,  2, 3, 4, 10, 14, 15);
	/* (a0 + a1 + a4 + a6 + a7 + a8)(b0 + b1 + b4 + b6 + b7 + b8)*(X^2 + X^4 + X0 + X1 + X2 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a4, a6, a7, a8),
				mzd_sum(t2, 6, b0, b1, b4, b6, b7, b8), 0), X, 6,  2, 4, 10, 11, 12, 15);
	/* (a0 + a2 + a3 + a4 + a7)(b0 + b2 + b3 + b4 + b7)*(X^2 + X0 + X1 + X3) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a3, a4, a7),
				mzd_sum(t2, 5, b0, b2, b3, b4, b7), 0), X, 4,  2, 10, 11, 13);
	/* (a1 + a2 + a3 + a6 + a8)(b1 + b2 + b3 + b6 + b8)*(X^3 + X1 + X2 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a2, a3, a6, a8),
				mzd_sum(t2, 5, b1, b2, b3, b6, b8), 0), X, 4,  3, 11, 12, 14);
	/* (a0 + a3 + a4 + a5 + a7)(b0 + b3 + b4 + b5 + b7)*(X^3 + X^4 + X^5 + X^6 + X^9 + X2 + X3 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a3, a4, a5, a7),
				mzd_sum(t2, 5, b0, b3, b4, b5, b7), 0), X, 8,  3, 4, 5, 6, 9, 12, 13, 14);
	/* (a1 + a4 + a5 + a6 + a8)(b1 + b4 + b5 + b6 + b8)*(X^4 + X^5 + X^6 + X^7 + X0 + X3 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a4, a5, a6, a8),
				mzd_sum(t2, 5, b1, b4, b5, b6, b8), 0), X, 8,  4, 5, 6, 7, 10, 13, 14, 15);
	/* (a2 + a3 + a4 + a6)(b2 + b3 + b4 + b6)*(X^2 + X^3 + X^6 + X^7 + X^8 + X0 + X1 + X2 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a2, a3, a4, a6),
				mzd_sum(t2, 4, b2, b3, b4, b6), 0), X, 10,  2, 3, 6, 7, 8, 10, 11, 12, 14, 15);
	/* (a0 + a1 + a3 + a6 + a7 + a8)(b0 + b1 + b3 + b6 + b7 + b8)*(X^2 + X^6 + X^8 + X^9 + X1 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a3, a6, a7, a8),
				mzd_sum(t2, 6, b0, b1, b3, b6, b7, b8), 0), X, 6,  2, 6, 8, 9, 11, 14);
	/* (a0 + a2 + a5 + a6 + a7)(b0 + b2 + b5 + b6 + b7)*(X^3 + X^7 + X^9 + X0 + X2 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a5, a6, a7),
				mzd_sum(t2, 5, b0, b2, b5, b6, b7), 0), X, 6,  3, 7, 9, 10, 12, 15);
	/* (a1 + a2 + a3 + a5 + a8)(b1 + b2 + b3 + b5 + b8)*(X^2 + X^3 + X^4 + X^5 + X^8 + X1 + X2 + X3) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a2, a3, a5, a8),
				mzd_sum(t2, 5, b1, b2, b3, b5, b8), 0), X, 8,  2, 3, 4, 5, 8, 11, 12, 13);
	/* (a0 + a4 + a7 + a8)(b0 + b4 + b7 + b8)*(X^4 + X^7 + X1 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a4, a7, a8),
				mzd_sum(t2, 4, b0, b4, b7, b8), 0), X, 4,  4, 7, 11, 14);
	/* (a1 + a4 + a5 + a7)(b1 + b4 + b5 + b7)*(X^5 + X^8 + X2 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a1, a4, a5, a7),
				mzd_sum(t2, 4, b1, b4, b5, b7), 0), X, 4,  5, 8, 12, 15);
	/* (a2 + a5 + a6 + a8)(b2 + b5 + b6 + b8)*(X^2 + X^3 + X^5 + X^6 + X^9 + X0 + X2 + X3) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a2, a5, a6, a8),
				mzd_sum(t2, 4, b2, b5, b6, b8), 0), X, 8,  2, 3, 5, 6, 9, 10, 12, 13);
	/* (a3 + a6 + a7)(b3 + b6 + b7)*(X^3 + X^4 + X^6 + X^7 + X0 + X1 + X3 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 3, a3, a6, a7),
				mzd_sum(t2, 3, b3, b6, b7), 0), X, 8,  3, 4, 6, 7, 10, 11, 13, 14);
	/* (a0 + a1 + a2 + a3 + a7)(b0 + b1 + b2 + b3 + b7)*(X^2 + X^4 + X^5 + X^7 + X^9 + X1 + X2 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a2, a3, a7),
				mzd_sum(t2, 5, b0, b1, b2, b3, b7), 0), X, 8,  2, 4, 5, 7, 9, 11, 12, 14);
	/* (a0 + a1 + a5 + a8)(b0 + b1 + b5 + b8)*(X^2 + X^8 + X^9 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a5, a8),
				mzd_sum(t2, 4, b0, b1, b5, b8), 0), X, 4,  2, 8, 9, 15);
	/* (a0 + a2 + a4 + a5 + a6 + a7)(b0 + b2 + b4 + b5 + b6 + b7)*(X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X0 + X1 + X2 + X3 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a4, a5, a6, a7),
				mzd_sum(t2, 6, b0, b2, b4, b5, b6, b7), 0), X, 12,  3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15);
	/* (a1 + a3 + a4 + a5 + a6)(b1 + b3 + b4 + b5 + b6)*(X^2 + X^3 + X^4 + X^6 + X^7 + X^8 + X^9 + X0 + X1 + X3 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a3, a4, a5, a6),
				mzd_sum(t2, 5, b1, b3, b4, b5, b6), 0), X, 12,  2, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15);
	/* (a2 + a3 + a5 + a7 + a8)(b2 + b3 + b5 + b7 + b8)*(X^3 + X^6 + X0 + X3) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a3, a5, a7, a8),
				mzd_sum(t2, 5, b2, b3, b5, b7, b8), 0), X, 4,  3, 6, 10, 13);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba10(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 36 multiplications and three temporaries (not optimal!) */
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];
	const mzd_t *a4 = A[4];
	const mzd_t *a5 = A[5];
	const mzd_t *a6 = A[6];
	const mzd_t *a7 = A[7];
	const mzd_t *a8 = A[8];
	const mzd_t *a9 = A[9];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];
	const mzd_t *b4 = B[4];
	const mzd_t *b5 = B[5];
	const mzd_t *b6 = B[6];
	const mzd_t *b7 = B[7];
	const mzd_t *b8 = B[8];
	const mzd_t *b9 = B[9];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	/* (a9)(b9)*(X + X^6 + X^8 + X1 + X3 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a9, b9, 0), X, 6,  1, 6, 8, 11, 13, 18);
	/* (a0)(b0)*(1 + X^5 + X^7 + X0 + X2 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a0, b0, 0), X, 6,  0, 5, 7, 10, 12, 17);
	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9)(b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9)*(X + X^2 + X^3 + X^4 + X^5 + X^8 + X^9 + X0 + X3 + X4 + X5 + X6 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9),
				mzd_sum(t2, 10, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9), 0), X, 13,  1, 2, 3, 4, 5, 8, 9, 10, 13, 14, 15, 16, 17);
	/* (a0 + a2 + a3 + a5 + a6 + a8 + a9)(b0 + b2 + b3 + b5 + b6 + b8 + b9)*(X^2 + X^3 + X^5 + X^6 + X^7 + X2 + X3 + X4 + X6 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a3, a5, a6, a8, a9),
				mzd_sum(t2, 7, b0, b2, b3, b5, b6, b8, b9), 0), X, 10,  2, 3, 5, 6, 7, 12, 13, 14, 16, 17);
	/* (a1 + a2 + a4 + a5 + a7 + a8)(b1 + b2 + b4 + b5 + b7 + b8)*(X + X^3 + X^4 + X^7 + X1 + X4 + X5 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a2, a4, a5, a7, a8),
				mzd_sum(t2, 6, b1, b2, b4, b5, b7, b8), 0), X, 8,  1, 3, 4, 7, 11, 14, 15, 17);
	/* (a0 + a1 + a3 + a4 + a6 + a7 + a9)(b0 + b1 + b3 + b4 + b6 + b7 + b9)*(X + X^2 + X^4 + X^5 + X^6 + X1 + X2 + X3 + X5 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a3, a4, a6, a7, a9),
				mzd_sum(t2, 7, b0, b1, b3, b4, b6, b7, b9), 0), X, 10,  1, 2, 4, 5, 6, 11, 12, 13, 15, 16);
	/* (a0 + a3 + a5 + a6 + a7)(b0 + b3 + b5 + b6 + b7)*(X + X^3 + X^7 + X^8 + X2 + X3 + X6 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a3, a5, a6, a7),
				mzd_sum(t2, 5, b0, b3, b5, b6, b7), 0), X, 8,  1, 3, 7, 8, 12, 13, 16, 17);
	/* (a1 + a3 + a4 + a5 + a8)(b1 + b3 + b4 + b5 + b8)*(X + X^2 + X^4 + X^6 + X^9 + X1 + X4 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a3, a4, a5, a8),
				mzd_sum(t2, 5, b1, b3, b4, b5, b8), 0), X, 8,  1, 2, 4, 6, 9, 11, 14, 17);
	/* (a2 + a4 + a5 + a6 + a9)(b2 + b4 + b5 + b6 + b9)*(X + X^2 + X^3 + X^5 + X^6 + X^7 + X^8 + X0 + X1 + X2 + X3 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a4, a5, a6, a9),
				mzd_sum(t2, 5, b2, b4, b5, b6, b9), 0), X, 12,  1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15);
	/* (a0 + a1 + a4 + a6 + a7 + a8)(b0 + b1 + b4 + b6 + b7 + b8)*(X^2 + X^3 + X^4 + X^6 + X^7 + X^8 + X^9 + X1 + X2 + X3 + X4 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a4, a6, a7, a8),
				mzd_sum(t2, 6, b0, b1, b4, b6, b7, b8), 0), X, 12,  2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16);
	/* (a0 + a2 + a3 + a4 + a7 + a9)(b0 + b2 + b3 + b4 + b7 + b9)*(X^3 + X^4 + X^5 + X^7 + X^8 + X^9 + X0 + X2 + X3 + X4 + X5 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a3, a4, a7, a9),
				mzd_sum(t2, 6, b0, b2, b3, b4, b7, b9), 0), X, 12,  3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 15, 17);
	/* (a1 + a2 + a3 + a6 + a8 + a9)(b1 + b2 + b3 + b6 + b8 + b9)*(X + X^4 + X^5 + X^9 + X0 + X4 + X5 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a2, a3, a6, a8, a9),
				mzd_sum(t2, 6, b1, b2, b3, b6, b8, b9), 0), X, 8,  1, 4, 5, 9, 10, 14, 15, 16);
	/* (a0 + a3 + a4 + a5 + a7)(b0 + b3 + b4 + b5 + b7)*(X^3 + X^5 + X^6 + X^7 + X^8 + X0 + X1 + X2 + X3 + X5 + X6 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a3, a4, a5, a7),
				mzd_sum(t2, 5, b0, b3, b4, b5, b7), 0), X, 12,  3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17);
	/* (a1 + a4 + a5 + a6 + a8)(b1 + b4 + b5 + b6 + b8)*(X + X^4 + X^7 + X^9 + X2 + X4 + X6 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a4, a5, a6, a8),
				mzd_sum(t2, 5, b1, b4, b5, b6, b8), 0), X, 8,  1, 4, 7, 9, 12, 14, 16, 17);
	/* (a2 + a3 + a4 + a6 + a9)(b2 + b3 + b4 + b6 + b9)*(X + X^2 + X^5 + X^6 + X0 + X1 + X5 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a3, a4, a6, a9),
				mzd_sum(t2, 5, b2, b3, b4, b6, b9), 0), X, 8,  1, 2, 5, 6, 10, 11, 15, 17);
	/* (a0 + a1 + a3 + a6 + a7 + a8)(b0 + b1 + b3 + b6 + b7 + b8)*(X^2 + X^3 + X^4 + X^8 + X^9 + X3 + X4 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a3, a6, a7, a8),
				mzd_sum(t2, 6, b0, b1, b3, b6, b7, b8), 0), X, 8,  2, 3, 4, 8, 9, 13, 14, 17);
	/* (a0 + a2 + a5 + a6 + a7 + a9)(b0 + b2 + b5 + b6 + b7 + b9)*(X + X^3 + X^4 + X^5 + X^6 + X^8 + X^9 + X0 + X1 + X3 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a5, a6, a7, a9),
				mzd_sum(t2, 6, b0, b2, b5, b6, b7, b9), 0), X, 12,  1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15);
	/* (a1 + a2 + a3 + a5 + a8 + a9)(b1 + b2 + b3 + b5 + b8 + b9)*(X^2 + X^4 + X^5 + X^6 + X^7 + X^9 + X0 + X1 + X2 + X4 + X5 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a2, a3, a5, a8, a9),
				mzd_sum(t2, 6, b1, b2, b3, b5, b8, b9), 0), X, 12,  2, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16);
	/* (a0 + a4 + a7 + a8)(b0 + b4 + b7 + b8)*(X + X^2 + X^4 + X^8 + X0 + X1 + X6 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a4, a7, a8),
				mzd_sum(t2, 4, b0, b4, b7, b8), 0), X, 8,  1, 2, 4, 8, 10, 11, 16, 17);
	/* (a1 + a4 + a5 + a7 + a9)(b1 + b4 + b5 + b7 + b9)*(X + X^2 + X^3 + X^5 + X^6 + X^8 + X^9 + X2 + X3 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a4, a5, a7, a9),
				mzd_sum(t2, 5, b1, b4, b5, b7, b9), 0), X, 10,  1, 2, 3, 5, 6, 8, 9, 12, 13, 17);
	/* (a2 + a5 + a6 + a8)(b2 + b5 + b6 + b8)*(X + X^2 + X^3 + X^4 + X^7 + X^8 + X^9 + X0 + X1 + X4) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a2, a5, a6, a8),
				mzd_sum(t2, 4, b2, b5, b6, b8), 0), X, 10,  1, 2, 3, 4, 7, 8, 9, 10, 11, 14);
	/* (a3 + a6 + a7 + a9)(b3 + b6 + b7 + b9)*(X^2 + X^3 + X^4 + X^5 + X^8 + X^9 + X0 + X1 + X2 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a3, a6, a7, a9),
				mzd_sum(t2, 4, b3, b6, b7, b9), 0), X, 10,  2, 3, 4, 5, 8, 9, 10, 11, 12, 15);
	/* (a0 + a1 + a2 + a3 + a7)(b0 + b1 + b2 + b3 + b7)*(X + X^5 + X^7 + X2 + X4 + X5) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a2, a3, a7),
				mzd_sum(t2, 5, b0, b1, b2, b3, b7), 0), X, 6,  1, 5, 7, 12, 14, 15);
	/* (a0 + a1 + a5 + a8 + a9)(b0 + b1 + b5 + b8 + b9)*(X + X^3 + X^4 + X^6 + X^7 + X^9 + X0 + X1 + X3 + X4 + X5 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a5, a8, a9),
				mzd_sum(t2, 5, b0, b1, b5, b8, b9), 0), X, 12,  1, 3, 4, 6, 7, 9, 10, 11, 13, 14, 15, 16);
	/* (a0 + a2 + a4 + a5 + a6 + a7)(b0 + b2 + b4 + b5 + b6 + b7)*(X + X^4 + X^6 + X0 + X1 + X3 + X5 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a4, a5, a6, a7),
				mzd_sum(t2, 6, b0, b2, b4, b5, b6, b7), 0), X, 8,  1, 4, 6, 10, 11, 13, 15, 17);
	/* (a1 + a3 + a4 + a5 + a6)(b1 + b3 + b4 + b5 + b6)*(X + X^2 + X^5 + X^6 + X^7 + X^8 + X2 + X3 + X4 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a3, a4, a5, a6),
				mzd_sum(t2, 5, b1, b3, b4, b5, b6), 0), X, 10,  1, 2, 5, 6, 7, 8, 12, 13, 14, 16);
	/* (a2 + a3 + a5 + a7 + a8 + a9)(b2 + b3 + b5 + b7 + b8 + b9)*(X + X^3 + X^5 + X^9 + X2 + X5 + X6 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a3, a5, a7, a8, a9),
				mzd_sum(t2, 6, b2, b3, b5, b7, b8, b9), 0), X, 8,  1, 3, 5, 9, 12, 15, 16, 17);
	/* (a0 + a4 + a5 + a6 + a7 + a9)(b0 + b4 + b5 + b6 + b7 + b9)*(X + X^2 + X^7 + X^8 + X0 + X4 + X6 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a4, a5, a6, a7, a9),
				mzd_sum(t2, 6, b0, b4, b5, b6, b7, b9), 0), X, 8,  1, 2, 7, 8, 10, 14, 16, 17);
	/* (a1 + a5 + a6 + a7 + a8)(b1 + b5 + b6 + b7 + b8)*(X + X^2 + X^3 + X^6 + X^9 + X3 + X5 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a5, a6, a7, a8),
				mzd_sum(t2, 5, b1, b5, b6, b7, b8), 0), X, 8,  1, 2, 3, 6, 9, 13, 15, 17);
	/* (a2 + a6 + a7 + a8 + a9)(b2 + b6 + b7 + b8 + b9)*(X + X^2 + X^3 + X^4 + X^6 + X^7 + X^8 + X0 + X1 + X3 + X4 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a6, a7, a8, a9),
				mzd_sum(t2, 5, b2, b6, b7, b8, b9), 0), X, 12,  1, 2, 3, 4, 6, 7, 8, 10, 11, 13, 14, 16);
	/* (a3 + a4 + a5 + a6 + a8)(b3 + b4 + b5 + b6 + b8)*(X^2 + X^3 + X^4 + X^5 + X^7 + X^8 + X^9 + X1 + X2 + X4 + X5 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a3, a4, a5, a6, a8),
				mzd_sum(t2, 5, b3, b4, b5, b6, b8), 0), X, 12,  2, 3, 4, 5, 7, 8, 9, 11, 12, 14, 15, 17);
	/* (a0 + a1 + a2 + a3 + a5 + a7 + a8)(b0 + b1 + b2 + b3 + b5 + b7 + b8)*(X^3 + X^6 + X^7 + X^8 + X^9 + X0 + X3 + X4 + X5 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a2, a3, a5, a7, a8),
				mzd_sum(t2, 7, b0, b1, b2, b3, b5, b7, b8), 0), X, 10,  3, 6, 7, 8, 9, 10, 13, 14, 15, 16);
	/* (a0 + a1 + a4 + a8 + a9)(b0 + b1 + b4 + b8 + b9)*(X + X^3 + X^4 + X^5 + X^9 + X0 + X1 + X2 + X5 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a4, a8, a9),
				mzd_sum(t2, 5, b0, b1, b4, b8, b9), 0), X, 10,  1, 3, 4, 5, 9, 10, 11, 12, 15, 16);
	/* (a0 + a2 + a4 + a5 + a8)(b0 + b2 + b4 + b5 + b8)*(X^2 + X^3 + X^5 + X0 + X2 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a4, a5, a8),
				mzd_sum(t2, 5, b0, b2, b4, b5, b8), 0), X, 6,  2, 3, 5, 10, 12, 16);
	/* (a1 + a3 + a4 + a7)(b1 + b3 + b4 + b7)*(X^3 + X^4 + X^6 + X1 + X3 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a1, a3, a4, a7),
				mzd_sum(t2, 4, b1, b3, b4, b7), 0), X, 6,  3, 4, 6, 11, 13, 17);
	/* (a2 + a3 + a4 + a5 + a7 + a9)(b2 + b3 + b4 + b5 + b7 + b9)*(X + X^3 + X^5 + X^7 + X^8 + X2 + X4 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a3, a4, a5, a7, a9),
				mzd_sum(t2, 6, b2, b3, b4, b5, b7, b9), 0), X, 8,  1, 3, 5, 7, 8, 12, 14, 17);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba11(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 40 multiplications and three temporaries (not optimal!) */
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];
	const mzd_t *a4 = A[4];
	const mzd_t *a5 = A[5];
	const mzd_t *a6 = A[6];
	const mzd_t *a7 = A[7];
	const mzd_t *a8 = A[8];
	const mzd_t *a9 = A[9];
	const mzd_t *a10 = A[10];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];
	const mzd_t *b4 = B[4];
	const mzd_t *b5 = B[5];
	const mzd_t *b6 = B[6];
	const mzd_t *b7 = B[7];
	const mzd_t *b8 = B[8];
	const mzd_t *b9 = B[9];
	const mzd_t *b10 = B[10];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	/* (a10)(b10)*(X^2 + X^3 + X^7 + X^8 + X^9 + X0 + X2 + X3 + X4 + X5 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a10, b10, 0), X, 12,  2, 3, 7, 8, 9, 10, 12, 13, 14, 15, 19, 20);
	/* (a0 + a2 + a3 + a5 + a6 + a8 + a9)(b0 + b2 + b3 + b5 + b6 + b8 + b9)*(X^3 + X^5 + X^6 + X^9 + X3 + X6 + X7 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a3, a5, a6, a8, a9),
				mzd_sum(t2, 7, b0, b2, b3, b5, b6, b8, b9), 0), X, 8,  3, 5, 6, 9, 13, 16, 17, 19);
	/* (a1 + a2 + a4 + a5 + a7 + a8 + a10)(b1 + b2 + b4 + b5 + b7 + b8 + b10)*(X^2 + X^3 + X^4 + X^6 + X^8 + X^9 + X2 + X3 + X5 + X7 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a2, a4, a5, a7, a8, a10),
				mzd_sum(t2, 7, b1, b2, b4, b5, b7, b8, b10), 0), X, 12,  2, 3, 4, 6, 8, 9, 12, 13, 15, 17, 18, 19);
	/* (a0 + a1 + a3 + a4 + a6 + a7 + a9 + a10)(b0 + b1 + b3 + b4 + b6 + b7 + b9 + b10)*(X^2 + X^4 + X^5 + X^8 + X2 + X5 + X6 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a3, a4, a6, a7, a9, a10),
				mzd_sum(t2, 8, b0, b1, b3, b4, b6, b7, b9, b10), 0), X, 8,  2, 4, 5, 8, 12, 15, 16, 18);
	/* (a0)(b0)*(1 + X + X^5 + X^6 + X^7 + X^8 + X0 + X1 + X2 + X3 + X7 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a0, b0, 0), X, 12,  0, 1, 5, 6, 7, 8, 10, 11, 12, 13, 17, 18);
	/* (a1)(b1)*(X + X^2 + X^6 + X^7 + X^8 + X^9 + X1 + X2 + X3 + X4 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a1, b1, 0), X, 12,  1, 2, 6, 7, 8, 9, 11, 12, 13, 14, 18, 19);
	/* (a0 + a1)(b0 + b1)*(X + X^2 + X^6 + X^7 + X^8 + X^9 + X1 + X2 + X3 + X4 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 12,  1, 2, 6, 7, 8, 9, 11, 12, 13, 14, 18, 19);
	/* (a0 + a2 + a4 + a6 + a8 + a10)(b0 + b2 + b4 + b6 + b8 + b10)*(X^2 + X^7 + X^9 + X2 + X4 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a4, a6, a8, a10),
				mzd_sum(t2, 6, b0, b2, b4, b6, b8, b10), 0), X, 6,  2, 7, 9, 12, 14, 19);
	/* (a1 + a3 + a5 + a7 + a9)(b1 + b3 + b5 + b7 + b9)*(X^2 + X^7 + X^9 + X2 + X4 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a3, a5, a7, a9),
				mzd_sum(t2, 5, b1, b3, b5, b7, b9), 0), X, 6,  2, 7, 9, 12, 14, 19);
	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10)(b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10)*(X^3 + X^4 + X^5 + X^6 + X^7 + X0 + X1 + X2 + X5 + X6 + X7 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 11, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10),
				mzd_sum(t2, 11, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10), 0), X, 13,  3, 4, 5, 6, 7, 10, 11, 12, 15, 16, 17, 18, 19);
	/* (a0 + a3 + a5 + a6 + a7 + a10)(b0 + b3 + b5 + b6 + b7 + b10)*(X^3 + X^6 + X^7 + X1 + X2 + X6 + X7 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a3, a5, a6, a7, a10),
				mzd_sum(t2, 6, b0, b3, b5, b6, b7, b10), 0), X, 8,  3, 6, 7, 11, 12, 16, 17, 18);
	/* (a1 + a3 + a4 + a5 + a8 + a10)(b1 + b3 + b4 + b5 + b8 + b10)*(X^4 + X^7 + X^8 + X2 + X3 + X7 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a4, a5, a8, a10),
				mzd_sum(t2, 6, b1, b3, b4, b5, b8, b10), 0), X, 8,  4, 7, 8, 12, 13, 17, 18, 19);
	/* (a2 + a4 + a5 + a6 + a9)(b2 + b4 + b5 + b6 + b9)*(X^2 + X^3 + X^5 + X^7 + X0 + X2 + X5 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a4, a5, a6, a9),
				mzd_sum(t2, 5, b2, b4, b5, b6, b9), 0), X, 8,  2, 3, 5, 7, 10, 12, 15, 18);
	/* (a0 + a1 + a4 + a6 + a7 + a8)(b0 + b1 + b4 + b6 + b7 + b8)*(X^3 + X^4 + X^6 + X^8 + X1 + X3 + X6 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a4, a6, a7, a8),
				mzd_sum(t2, 6, b0, b1, b4, b6, b7, b8), 0), X, 8,  3, 4, 6, 8, 11, 13, 16, 19);
	/* (a0 + a2 + a3 + a4 + a7 + a9 + a10)(b0 + b2 + b3 + b4 + b7 + b9 + b10)*(X^2 + X^3 + X^4 + X^5 + X^8 + X0 + X3 + X5 + X7 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a3, a4, a7, a9, a10),
				mzd_sum(t2, 7, b0, b2, b3, b4, b7, b9, b10), 0), X, 10,  2, 3, 4, 5, 8, 10, 13, 15, 17, 19);
	/* (a1 + a2 + a3 + a6 + a8 + a9 + a10)(b1 + b2 + b3 + b6 + b8 + b9 + b10)*(X^2 + X^4 + X^5 + X^6 + X^7 + X^8 + X0 + X1 + X2 + X3 + X5 + X6 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a2, a3, a6, a8, a9, a10),
				mzd_sum(t2, 7, b1, b2, b3, b6, b8, b9, b10), 0), X, 14,  2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 18, 19);
	/* (a0 + a3 + a4 + a5 + a7 + a10)(b0 + b3 + b4 + b5 + b7 + b10)*(X^2 + X^3 + X^5 + X^6 + X^8 + X^9 + X0 + X1 + X3 + X4 + X5 + X6 + X7 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a3, a4, a5, a7, a10),
				mzd_sum(t2, 6, b0, b3, b4, b5, b7, b10), 0), X, 14,  2, 3, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19);
	/* (a1 + a4 + a5 + a6 + a8)(b1 + b4 + b5 + b6 + b8)*(X^2 + X^4 + X^6 + X^8 + X1 + X3 + X6 + X7 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a4, a5, a6, a8),
				mzd_sum(t2, 5, b1, b4, b5, b6, b8), 0), X, 10,  2, 4, 6, 8, 11, 13, 16, 17, 18, 19);
	/* (a2 + a3 + a4 + a6 + a9 + a10)(b2 + b3 + b4 + b6 + b9 + b10)*(X^2 + X^5 + X^8 + X0 + X3 + X5 + X7 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a3, a4, a6, a9, a10),
				mzd_sum(t2, 6, b2, b3, b4, b6, b9, b10), 0), X, 8,  2, 5, 8, 10, 13, 15, 17, 18);
	/* (a0 + a1 + a3 + a6 + a7 + a8 + a10)(b0 + b1 + b3 + b6 + b7 + b8 + b10)*(X^2 + X^3 + X^4 + X^8 + X^9 + X3 + X4 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a3, a6, a7, a8, a10),
				mzd_sum(t2, 7, b0, b1, b3, b6, b7, b8, b10), 0), X, 8,  2, 3, 4, 8, 9, 13, 14, 17);
	/* (a0 + a2 + a5 + a6 + a7 + a9)(b0 + b2 + b5 + b6 + b7 + b9)*(X^3 + X^4 + X^5 + X^9 + X0 + X4 + X5 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a5, a6, a7, a9),
				mzd_sum(t2, 6, b0, b2, b5, b6, b7, b9), 0), X, 8,  3, 4, 5, 9, 10, 14, 15, 18);
	/* (a1 + a2 + a3 + a5 + a8 + a9 + a10)(b1 + b2 + b3 + b5 + b8 + b9 + b10)*(X^4 + X^5 + X^6 + X0 + X1 + X5 + X6 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a2, a3, a5, a8, a9, a10),
				mzd_sum(t2, 7, b1, b2, b3, b5, b8, b9, b10), 0), X, 8,  4, 5, 6, 10, 11, 15, 16, 19);
	/* (a0 + a4 + a7 + a8 + a10)(b0 + b4 + b7 + b8 + b10)*(X^2 + X^4 + X^6 + X0 + X3 + X6 + X7 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a4, a7, a8, a10),
				mzd_sum(t2, 5, b0, b4, b7, b8, b10), 0), X, 8,  2, 4, 6, 10, 13, 16, 17, 18);
	/* (a1 + a4 + a5 + a7 + a9 + a10)(b1 + b4 + b5 + b7 + b9 + b10)*(X^3 + X^5 + X^7 + X1 + X4 + X7 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a4, a5, a7, a9, a10),
				mzd_sum(t2, 6, b1, b4, b5, b7, b9, b10), 0), X, 8,  3, 5, 7, 11, 14, 17, 18, 19);
	/* (a2 + a5 + a6 + a8 + a10)(b2 + b5 + b6 + b8 + b10)*(X^2 + X^3 + X^4 + X^6 + X^7 + X^9 + X0 + X3 + X4 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a5, a6, a8, a10),
				mzd_sum(t2, 5, b2, b5, b6, b8, b10), 0), X, 10,  2, 3, 4, 6, 7, 9, 10, 13, 14, 18);
	/* (a3 + a6 + a7 + a9)(b3 + b6 + b7 + b9)*(X^3 + X^4 + X^5 + X^7 + X^8 + X0 + X1 + X4 + X5 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a3, a6, a7, a9),
				mzd_sum(t2, 4, b3, b6, b7, b9), 0), X, 10,  3, 4, 5, 7, 8, 10, 11, 14, 15, 19);
	/* (a0 + a1 + a2 + a3 + a7 + a10)(b0 + b1 + b2 + b3 + b7 + b10)*(X^2 + X^5 + X^6 + X^8 + X^9 + X1 + X3 + X5 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a2, a3, a7, a10),
				mzd_sum(t2, 6, b0, b1, b2, b3, b7, b10), 0), X, 10,  2, 5, 6, 8, 9, 11, 13, 15, 18, 19);
	/* (a0 + a1 + a5 + a8 + a9)(b0 + b1 + b5 + b8 + b9)*(X^3 + X^4 + X^7 + X^8 + X^9 + X0 + X4 + X5 + X6 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a5, a8, a9),
				mzd_sum(t2, 5, b0, b1, b5, b8, b9), 0), X, 10,  3, 4, 7, 8, 9, 10, 14, 15, 16, 18);
	/* (a0 + a2 + a4 + a5 + a6 + a7)(b0 + b2 + b4 + b5 + b6 + b7)*(X^4 + X^8 + X0 + X5 + X7 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a4, a5, a6, a7),
				mzd_sum(t2, 6, b0, b2, b4, b5, b6, b7), 0), X, 6,  4, 8, 10, 15, 17, 18);
	/* (a1 + a3 + a4 + a5 + a6 + a10)(b1 + b3 + b4 + b5 + b6 + b10)*(X^5 + X^9 + X1 + X6 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a4, a5, a6, a10),
				mzd_sum(t2, 6, b1, b3, b4, b5, b6, b10), 0), X, 6,  5, 9, 11, 16, 18, 19);
	/* (a2 + a3 + a5 + a7 + a8 + a9 + a10)(b2 + b3 + b5 + b7 + b8 + b9 + b10)*(X^2 + X^3 + X^5 + X^6 + X^7 + X^8 + X1 + X3 + X4 + X5 + X6 + X7 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a5, a7, a8, a9, a10),
				mzd_sum(t2, 7, b2, b3, b5, b7, b8, b9, b10), 0), X, 14,  2, 3, 5, 6, 7, 8, 11, 13, 14, 15, 16, 17, 18, 19);
	/* (a0 + a4 + a5 + a6 + a7 + a9)(b0 + b4 + b5 + b6 + b7 + b9)*(X^2 + X^6 + X^7 + X0 + X1 + X3 + X4 + X6 + X7 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a4, a5, a6, a7, a9),
				mzd_sum(t2, 6, b0, b4, b5, b6, b7, b9), 0), X, 10,  2, 6, 7, 10, 11, 13, 14, 16, 17, 18);
	/* (a1 + a5 + a6 + a7 + a8 + a10)(b1 + b5 + b6 + b7 + b8 + b10)*(X^3 + X^7 + X^8 + X1 + X2 + X4 + X5 + X7 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a5, a6, a7, a8, a10),
				mzd_sum(t2, 6, b1, b5, b6, b7, b8, b10), 0), X, 10,  3, 7, 8, 11, 12, 14, 15, 17, 18, 19);
	/* (a2 + a6 + a7 + a8 + a9)(b2 + b6 + b7 + b8 + b9)*(X^2 + X^3 + X^4 + X^7 + X0 + X4 + X6 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a6, a7, a8, a9),
				mzd_sum(t2, 5, b2, b6, b7, b8, b9), 0), X, 8,  2, 3, 4, 7, 10, 14, 16, 18);
	/* (a3 + a4 + a5 + a6 + a8 + a10)(b3 + b4 + b5 + b6 + b8 + b10)*(X^3 + X^4 + X^5 + X^8 + X1 + X5 + X7 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a3, a4, a5, a6, a8, a10),
				mzd_sum(t2, 6, b3, b4, b5, b6, b8, b10), 0), X, 8,  3, 4, 5, 8, 11, 15, 17, 19);
	/* (a0 + a1 + a2 + a3 + a5 + a7 + a8)(b0 + b1 + b2 + b3 + b5 + b7 + b8)*(X^2 + X^3 + X^6 + X^8 + X0 + X2 + X3 + X5 + X6 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a2, a3, a5, a7, a8),
				mzd_sum(t2, 7, b0, b1, b2, b3, b5, b7, b8), 0), X, 10,  2, 3, 6, 8, 10, 12, 13, 15, 16, 19);
	/* (a0 + a1 + a4 + a8 + a9 + a10)(b0 + b1 + b4 + b8 + b9 + b10)*(X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X0 + X3 + X4 + X5 + X6 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a4, a8, a9, a10),
				mzd_sum(t2, 6, b0, b1, b4, b8, b9, b10), 0), X, 14,  2, 3, 4, 5, 6, 7, 8, 10, 13, 14, 15, 16, 18, 19);
	/* (a0 + a2 + a4 + a5 + a8)(b0 + b2 + b4 + b5 + b8)*(X^2 + X^3 + X^5 + X0 + X2 + X6) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a4, a5, a8),
				mzd_sum(t2, 5, b0, b2, b4, b5, b8), 0), X, 6,  2, 3, 5, 10, 12, 16);
	/* (a1 + a3 + a4 + a7)(b1 + b3 + b4 + b7)*(X^3 + X^4 + X^6 + X1 + X3 + X7) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a1, a3, a4, a7),
				mzd_sum(t2, 4, b1, b3, b4, b7), 0), X, 6,  3, 4, 6, 11, 13, 17);
	/* (a2 + a3 + a4 + a5 + a7 + a9 + a10)(b2 + b3 + b4 + b5 + b7 + b9 + b10)*(X^3 + X^5 + X^6 + X^7 + X1 + X2 + X3 + X4 + X7 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a4, a5, a7, a9, a10),
				mzd_sum(t2, 7, b2, b3, b4, b5, b7, b9, b10), 0), X, 10,  3, 5, 6, 7, 11, 12, 13, 14, 17, 18);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}


void _mzd_ptr_addmul_karatsuba12(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 45 multiplications and three temporaries (not optimal!) */
	const mzd_t *a0 = A[0];
	const mzd_t *a1 = A[1];
	const mzd_t *a2 = A[2];
	const mzd_t *a3 = A[3];
	const mzd_t *a4 = A[4];
	const mzd_t *a5 = A[5];
	const mzd_t *a6 = A[6];
	const mzd_t *a7 = A[7];
	const mzd_t *a8 = A[8];
	const mzd_t *a9 = A[9];
	const mzd_t *a10 = A[10];
	const mzd_t *a11 = A[11];

	const mzd_t *b0 = B[0];
	const mzd_t *b1 = B[1];
	const mzd_t *b2 = B[2];
	const mzd_t *b3 = B[3];
	const mzd_t *b4 = B[4];
	const mzd_t *b5 = B[5];
	const mzd_t *b6 = B[6];
	const mzd_t *b7 = B[7];
	const mzd_t *b8 = B[8];
	const mzd_t *b9 = B[9];
	const mzd_t *b10 = B[10];
	const mzd_t *b11 = B[11];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);


	/* (a11)(b11)*(X + X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X6 + X7 + X8 + X9 + X^20 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a11, b11, 0), X, 14,  1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22);
	/* (a0)(b0)*(1 + X + X^2 + X^3 + X^4 + X^5 + X^6 + X5 + X6 + X7 + X8 + X9 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a0, b0, 0), X, 14,  0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21);
	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11)(b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10 + b11)*(X + X^3 + X^5 + X^7 + X^8 + X^9 + X0 + X1 + X2 + X3 + X4 + X5 + X7 + X9 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 12, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11),
				mzd_sum(t2, 12, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11), 0), X, 15,  1, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 21);
	/* (a0 + a2 + a3 + a5 + a6 + a8 + a9 + a11)(b0 + b2 + b3 + b5 + b6 + b8 + b9 + b11)*(X + X^2 + X^4 + X^5 + X^7 + X^9 + X0 + X2 + X3 + X5 + X7 + X8 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a3, a5, a6, a8, a9, a11),
				mzd_sum(t2, 8, b0, b2, b3, b5, b6, b8, b9, b11), 0), X, 14,  1, 2, 4, 5, 7, 9, 10, 12, 13, 15, 17, 18, 20, 21);
	/* (a1 + a2 + a4 + a5 + a7 + a8 + a10 + a11)(b1 + b2 + b4 + b5 + b7 + b8 + b10 + b11)*(X + X^4 + X^7 + X^8 + X0 + X1 + X3 + X4 + X7 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a2, a4, a5, a7, a8, a10, a11),
				mzd_sum(t2, 8, b1, b2, b4, b5, b7, b8, b10, b11), 0), X, 10,  1, 4, 7, 8, 10, 11, 13, 14, 17, 20);
	/* (a0 + a1 + a3 + a4 + a6 + a7 + a9 + a10)(b0 + b1 + b3 + b4 + b6 + b7 + b9 + b10)*(X^2 + X^5 + X^8 + X^9 + X1 + X2 + X4 + X5 + X8 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a3, a4, a6, a7, a9, a10),
				mzd_sum(t2, 8, b0, b1, b3, b4, b6, b7, b9, b10), 0), X, 10,  2, 5, 8, 9, 11, 12, 14, 15, 18, 21);
	/* (a0 + a3 + a5 + a6 + a7 + a10)(b0 + b3 + b5 + b6 + b7 + b10)*(X^2 + X^3 + X^4 + X^6 + X7 + X8 + X9 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a3, a5, a6, a7, a10),
				mzd_sum(t2, 6, b0, b3, b5, b6, b7, b10), 0), X, 8,  2, 3, 4, 6, 17, 18, 19, 21);
	/* (a1 + a3 + a4 + a5 + a8 + a10 + a11)(b1 + b3 + b4 + b5 + b8 + b10 + b11)*(X + X^2 + X^6 + X6 + X7 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a3, a4, a5, a8, a10, a11),
				mzd_sum(t2, 7, b1, b3, b4, b5, b8, b10, b11), 0), X, 6,  1, 2, 6, 16, 17, 21);
	/* (a2 + a4 + a5 + a6 + a9 + a11)(b2 + b4 + b5 + b6 + b9 + b11)*(X + X^4 + X^5 + X^6 + X6 + X9 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a4, a5, a6, a9, a11),
				mzd_sum(t2, 6, b2, b4, b5, b6, b9, b11), 0), X, 8,  1, 4, 5, 6, 16, 19, 20, 21);
	/* (a0 + a1 + a4 + a6 + a7 + a8 + a11)(b0 + b1 + b4 + b6 + b7 + b8 + b11)*(X + X^3 + X^4 + X6 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a4, a6, a7, a8, a11),
				mzd_sum(t2, 7, b0, b1, b4, b6, b7, b8, b11), 0), X, 6,  1, 3, 4, 16, 18, 19);
	/* (a0 + a2 + a3 + a4 + a7 + a9 + a10 + a11)(b0 + b2 + b3 + b4 + b7 + b9 + b10 + b11)*(X^2 + X^4 + X^5 + X7 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a3, a4, a7, a9, a10, a11),
				mzd_sum(t2, 8, b0, b2, b3, b4, b7, b9, b10, b11), 0), X, 6,  2, 4, 5, 17, 19, 20);
	/* (a1 + a2 + a3 + a6 + a8 + a9 + a10)(b1 + b2 + b3 + b6 + b8 + b9 + b10)*(X^3 + X^5 + X^6 + X8 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a2, a3, a6, a8, a9, a10),
				mzd_sum(t2, 7, b1, b2, b3, b6, b8, b9, b10), 0), X, 6,  3, 5, 6, 18, 20, 21);
	/* (a0 + a3 + a4 + a5 + a7 + a10 + a11)(b0 + b3 + b4 + b5 + b7 + b10 + b11)*(X^2 + X^4 + X^5 + X^6 + X7 + X9 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a3, a4, a5, a7, a10, a11),
				mzd_sum(t2, 7, b0, b3, b4, b5, b7, b10, b11), 0), X, 8,  2, 4, 5, 6, 17, 19, 20, 21);
	/* (a1 + a4 + a5 + a6 + a8 + a11)(b1 + b4 + b5 + b6 + b8 + b11)*(X + X^2 + X^4 + X6 + X7 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a4, a5, a6, a8, a11),
				mzd_sum(t2, 6, b1, b4, b5, b6, b8, b11), 0), X, 6,  1, 2, 4, 16, 17, 19);
	/* (a2 + a3 + a4 + a6 + a9 + a10 + a11)(b2 + b3 + b4 + b6 + b9 + b10 + b11)*(X^2 + X^3 + X^5 + X7 + X8 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a4, a6, a9, a10, a11),
				mzd_sum(t2, 7, b2, b3, b4, b6, b9, b10, b11), 0), X, 6,  2, 3, 5, 17, 18, 20);
	/* (a0 + a1 + a3 + a6 + a7 + a8 + a10)(b0 + b1 + b3 + b6 + b7 + b8 + b10)*(X + X^2 + X^3 + X^6 + X6 + X7 + X8 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a3, a6, a7, a8, a10),
				mzd_sum(t2, 7, b0, b1, b3, b6, b7, b8, b10), 0), X, 8,  1, 2, 3, 6, 16, 17, 18, 21);
	/* (a0 + a2 + a5 + a6 + a7 + a9)(b0 + b2 + b5 + b6 + b7 + b9)*(X + X^5 + X^6 + X6 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a5, a6, a7, a9),
				mzd_sum(t2, 6, b0, b2, b5, b6, b7, b9), 0), X, 6,  1, 5, 6, 16, 20, 21);
	/* (a1 + a2 + a3 + a5 + a8 + a9 + a10)(b1 + b2 + b3 + b5 + b8 + b9 + b10)*(X + X^3 + X^4 + X^5 + X6 + X8 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a2, a3, a5, a8, a9, a10),
				mzd_sum(t2, 7, b1, b2, b3, b5, b8, b9, b10), 0), X, 8,  1, 3, 4, 5, 16, 18, 19, 20);
	/* (a0 + a4 + a7 + a8 + a10)(b0 + b4 + b7 + b8 + b10)*(X + X^2 + X^3 + X^9 + X3 + X4 + X5 + X7 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a4, a7, a8, a10),
				mzd_sum(t2, 5, b0, b4, b7, b8, b10), 0), X, 10,  1, 2, 3, 9, 13, 14, 15, 17, 20, 21);
	/* (a1 + a4 + a5 + a7 + a9 + a10 + a11)(b1 + b4 + b5 + b7 + b9 + b10 + b11)*(X + X^5 + X^6 + X^7 + X0 + X4 + X5 + X7 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a4, a5, a7, a9, a10, a11),
				mzd_sum(t2, 7, b1, b4, b5, b7, b9, b10, b11), 0), X, 10,  1, 5, 6, 7, 10, 14, 15, 17, 19, 20);
	/* (a2 + a5 + a6 + a8 + a10 + a11)(b2 + b5 + b6 + b8 + b10 + b11)*(X^2 + X^6 + X^7 + X^8 + X1 + X5 + X6 + X8 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a5, a6, a8, a10, a11),
				mzd_sum(t2, 6, b2, b5, b6, b8, b10, b11), 0), X, 10,  2, 6, 7, 8, 11, 15, 16, 18, 20, 21);
	/* (a3 + a6 + a7 + a9 + a11)(b3 + b6 + b7 + b9 + b11)*(X + X^2 + X^4 + X^5 + X^6 + X^8 + X^9 + X2 + X8 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a3, a6, a7, a9, a11),
				mzd_sum(t2, 5, b3, b6, b7, b9, b11), 0), X, 10,  1, 2, 4, 5, 6, 8, 9, 12, 18, 20);
	/* (a0 + a1 + a2 + a3 + a7 + a10 + a11)(b0 + b1 + b2 + b3 + b7 + b10 + b11)*(X + X^4 + X^5 + X^7 + X^9 + X1 + X2 + X5 + X6 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a2, a3, a7, a10, a11),
				mzd_sum(t2, 7, b0, b1, b2, b3, b7, b10, b11), 0), X, 10,  1, 4, 5, 7, 9, 11, 12, 15, 16, 21);
	/* (a0 + a1 + a5 + a8 + a9 + a11)(b0 + b1 + b5 + b8 + b9 + b11)*(X + X^2 + X^3 + X^4 + X^6 + X0 + X1 + X2 + X3 + X5 + X6 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a5, a8, a9, a11),
				mzd_sum(t2, 6, b0, b1, b5, b8, b9, b11), 0), X, 12,  1, 2, 3, 4, 6, 10, 11, 12, 13, 15, 16, 19);
	/* (a0 + a2 + a4 + a5 + a6 + a7 + a11)(b0 + b2 + b4 + b5 + b6 + b7 + b11)*(X^2 + X^4 + X^7 + X^8 + X^9 + X0 + X2 + X4 + X5 + X7 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a4, a5, a6, a7, a11),
				mzd_sum(t2, 7, b0, b2, b4, b5, b6, b7, b11), 0), X, 12,  2, 4, 7, 8, 9, 10, 12, 14, 15, 17, 18, 19);
	/* (a1 + a3 + a4 + a5 + a6 + a10)(b1 + b3 + b4 + b5 + b6 + b10)*(X^3 + X^5 + X^8 + X^9 + X0 + X1 + X3 + X5 + X6 + X8 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a4, a5, a6, a10),
				mzd_sum(t2, 6, b1, b3, b4, b5, b6, b10), 0), X, 12,  3, 5, 8, 9, 10, 11, 13, 15, 16, 18, 19, 20);
	/* (a2 + a3 + a5 + a7 + a8 + a9 + a10)(b2 + b3 + b5 + b7 + b8 + b9 + b10)*(X^3 + X^4 + X^5 + X^6 + X^8 + X2 + X3 + X4 + X5 + X7 + X8 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a5, a7, a8, a9, a10),
				mzd_sum(t2, 7, b2, b3, b5, b7, b8, b9, b10), 0), X, 12,  3, 4, 5, 6, 8, 12, 13, 14, 15, 17, 18, 21);
	/* (a0 + a4 + a5 + a6 + a7 + a9 + a11)(b0 + b4 + b5 + b6 + b7 + b9 + b11)*(X + X^3 + X^9 + X2 + X3 + X5 + X6 + X7 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a4, a5, a6, a7, a9, a11),
				mzd_sum(t2, 7, b0, b4, b5, b6, b7, b9, b11), 0), X, 10,  1, 3, 9, 12, 13, 15, 16, 17, 19, 20);
	/* (a1 + a5 + a6 + a7 + a8 + a10)(b1 + b5 + b6 + b7 + b8 + b10)*(X^2 + X^4 + X0 + X3 + X4 + X6 + X7 + X8 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a5, a6, a7, a8, a10),
				mzd_sum(t2, 6, b1, b5, b6, b7, b8, b10), 0), X, 10,  2, 4, 10, 13, 14, 16, 17, 18, 20, 21);
	/* (a2 + a6 + a7 + a8 + a9 + a11)(b2 + b6 + b7 + b8 + b9 + b11)*(X + X^2 + X^4 + X^6 + X^7 + X1 + X4 + X5 + X6 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a6, a7, a8, a9, a11),
				mzd_sum(t2, 6, b2, b6, b7, b8, b9, b11), 0), X, 10,  1, 2, 4, 6, 7, 11, 14, 15, 16, 20);
	/* (a3 + a4 + a5 + a6 + a8 + a10 + a11)(b3 + b4 + b5 + b6 + b8 + b10 + b11)*(X^2 + X^3 + X^5 + X^7 + X^8 + X2 + X5 + X6 + X7 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a3, a4, a5, a6, a8, a10, a11),
				mzd_sum(t2, 7, b3, b4, b5, b6, b8, b10, b11), 0), X, 10,  2, 3, 5, 7, 8, 12, 15, 16, 17, 21);
	/* (a0 + a1 + a2 + a3 + a5 + a7 + a8 + a11)(b0 + b1 + b2 + b3 + b5 + b7 + b8 + b11)*(X + X^2 + X^3 + X^4 + X^9 + X0 + X2 + X4 + X5 + X8 + X9 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a2, a3, a5, a7, a8, a11),
				mzd_sum(t2, 8, b0, b1, b2, b3, b5, b7, b8, b11), 0), X, 12,  1, 2, 3, 4, 9, 10, 12, 14, 15, 18, 19, 21);
	/* (a0 + a1 + a4 + a8 + a9 + a10 + a11)(b0 + b1 + b4 + b8 + b9 + b10 + b11)*(X + X^2 + X^5 + X^7 + X^8 + X^9 + X3 + X9 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a4, a8, a9, a10, a11),
				mzd_sum(t2, 7, b0, b1, b4, b8, b9, b10, b11), 0), X, 10,  1, 2, 5, 7, 8, 9, 13, 19, 20, 21);
	/* (a0 + a2 + a4 + a5 + a8)(b0 + b2 + b4 + b5 + b8)*(X + X^2 + X^3 + X^5 + X^6 + X^8 + X0 + X1 + X2 + X3 + X6 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a4, a5, a8),
				mzd_sum(t2, 5, b0, b2, b4, b5, b8), 0), X, 12,  1, 2, 3, 5, 6, 8, 10, 11, 12, 13, 16, 18);
	/* (a1 + a3 + a4 + a7 + a11)(b1 + b3 + b4 + b7 + b11)*(X^2 + X^3 + X^4 + X^6 + X^7 + X^9 + X1 + X2 + X3 + X4 + X7 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a3, a4, a7, a11),
				mzd_sum(t2, 5, b1, b3, b4, b7, b11), 0), X, 12,  2, 3, 4, 6, 7, 9, 11, 12, 13, 14, 17, 19);
	/* (a2 + a3 + a4 + a5 + a7 + a9 + a10)(b2 + b3 + b4 + b5 + b7 + b9 + b10)*(X^2 + X^5 + X^6 + X^8 + X^9 + X0 + X1 + X5 + X7 + X8 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a4, a5, a7, a9, a10),
				mzd_sum(t2, 7, b2, b3, b4, b5, b7, b9, b10), 0), X, 12,  2, 5, 6, 8, 9, 10, 11, 15, 17, 18, 19, 20);
	/* (a0 + a4 + a5 + a9 + a10)(b0 + b4 + b5 + b9 + b10)*(X^2 + X^3 + X^4 + X^7 + X^8 + X0 + X1 + X2 + X3 + X5 + X6 + X9 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a4, a5, a9, a10),
				mzd_sum(t2, 5, b0, b4, b5, b9, b10), 0), X, 14,  2, 3, 4, 7, 8, 10, 11, 12, 13, 15, 16, 19, 20, 21);
	/* (a1 + a4 + a6 + a9 + a11)(b1 + b4 + b6 + b9 + b11)*(X + X^2 + X^6 + X^7 + X^8 + X^9 + X1 + X2 + X3 + X4 + X8 + X9) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a4, a6, a9, a11),
				mzd_sum(t2, 5, b1, b4, b6, b9, b11), 0), X, 12,  1, 2, 6, 7, 8, 9, 11, 12, 13, 14, 18, 19);
	/* (a2 + a4 + a7 + a9)(b2 + b4 + b7 + b9)*(X^2 + X^3 + X^7 + X^8 + X^9 + X0 + X2 + X3 + X4 + X5 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a2, a4, a7, a9),
				mzd_sum(t2, 4, b2, b4, b7, b9), 0), X, 12,  2, 3, 7, 8, 9, 10, 12, 13, 14, 15, 19, 20);
	/* (a3 + a4 + a8 + a9)(b3 + b4 + b8 + b9)*(X^3 + X^4 + X^8 + X^9 + X0 + X1 + X3 + X4 + X5 + X6 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a3, a4, a8, a9),
				mzd_sum(t2, 4, b3, b4, b8, b9), 0), X, 12,  3, 4, 8, 9, 10, 11, 13, 14, 15, 16, 20, 21);
	/* (a0 + a1 + a2 + a3 + a5 + a6 + a7 + a8 + a10 + a11)(b0 + b1 + b2 + b3 + b5 + b6 + b7 + b8 + b10 + b11)*(X + X^2 + X^3 + X^6 + X^7 + X^9 + X0 + X1 + X2 + X4 + X5 + X8 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a1, a2, a3, a5, a6, a7, a8, a10, a11),
				mzd_sum(t2, 10, b0, b1, b2, b3, b5, b6, b7, b8, b10, b11), 0), X, 14,  1, 2, 3, 6, 7, 9, 10, 11, 12, 14, 15, 18, 19, 20);
	/* (a0 + a1 + a5 + a6 + a10 + a11)(b0 + b1 + b5 + b6 + b10 + b11)*(X + X^6 + X^8 + X1 + X3 + X8) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a5, a6, a10, a11),
				mzd_sum(t2, 6, b0, b1, b5, b6, b10, b11), 0), X, 6,  1, 6, 8, 11, 13, 18);
	/* (a0 + a2 + a5 + a7 + a10)(b0 + b2 + b5 + b7 + b10)*(X + X^2 + X^4 + X^6 + X^7 + X^8 + X2 + X3 + X6 + X8 + X9 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a5, a7, a10),
				mzd_sum(t2, 5, b0, b2, b5, b7, b10), 0), X, 12,  1, 2, 4, 6, 7, 8, 12, 13, 16, 18, 19, 21);
	/* (a1 + a3 + a6 + a8 + a11)(b1 + b3 + b6 + b8 + b11)*(X + X^4 + X^6 + X^8 + X^9 + X3 + X4 + X6 + X8 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a3, a6, a8, a11),
				mzd_sum(t2, 5, b1, b3, b6, b8, b11), 0), X, 10,  1, 4, 6, 8, 9, 13, 14, 16, 18, 21);
	/* (a2 + a3 + a7 + a8)(b2 + b3 + b7 + b8)*(X^3 + X^8 + X0 + X3 + X5 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a2, a3, a7, a8),
				mzd_sum(t2, 4, b2, b3, b7, b8), 0), X, 6,  3, 8, 10, 13, 15, 20);
	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba13(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 49 multiplications and three temporaries (not optimal!) */
	const mzd_t *a0  = A[ 0];
	const mzd_t *a1  = A[ 1];
	const mzd_t *a2  = A[ 2];
	const mzd_t *a3  = A[ 3];
	const mzd_t *a4  = A[ 4];
	const mzd_t *a5  = A[ 5];
	const mzd_t *a6  = A[ 6];
	const mzd_t *a7  = A[ 7];
	const mzd_t *a8  = A[ 8];
	const mzd_t *a9  = A[ 9];
	const mzd_t *a10 = A[10];
	const mzd_t *a11 = A[11];
	const mzd_t *a12 = A[12];

	const mzd_t *b0  = B[ 0];
	const mzd_t *b1  = B[ 1];
	const mzd_t *b2  = B[ 2];
	const mzd_t *b3  = B[ 3];
	const mzd_t *b4  = B[ 4];
	const mzd_t *b5  = B[ 5];
	const mzd_t *b6  = B[ 6];
	const mzd_t *b7  = B[ 7];
	const mzd_t *b8  = B[ 8];
	const mzd_t *b9  = B[ 9];
	const mzd_t *b10 = B[10];
	const mzd_t *b11 = B[11];
	const mzd_t *b12 = B[12];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	/* (a12)(b12)*(X^2 + X^9 + X7 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a12, b12, 0), X, 4,  2, 9, 17, 24);
	/* (a0 + a2 + a3 + a5 + a6 + a8 + a9 + a11 + a12)(b0 + b2 + b3 + b5 + b6 + b8 + b9 + b11 + b12)*(X^2 + X^4 + X^5 + X^7 + X^8 + X^9 + X0 + X2 + X3 + X5 + X6 + X7 + X8 + X^20 + X^21 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a3, a5, a6, a8, a9, a11, a12),
				mzd_sum(t2, 9, b0, b2, b3, b5, b6, b8, b9, b11, b12), 0), X, 16,  2, 4, 5, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 20, 21, 23);
	/* (a1 + a2 + a4 + a5 + a7 + a8 + a10 + a11)(b1 + b2 + b4 + b5 + b7 + b8 + b10 + b11)*(X^2 + X^3 + X^5 + X^6 + X^8 + X0 + X1 + X3 + X4 + X6 + X8 + X9 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a2, a4, a5, a7, a8, a10, a11),
				mzd_sum(t2, 8, b1, b2, b4, b5, b7, b8, b10, b11), 0), X, 14,  2, 3, 5, 6, 8, 10, 11, 13, 14, 16, 18, 19, 21, 22);
	/* (a0 + a1 + a3 + a4 + a6 + a7 + a9 + a10 + a12)(b0 + b1 + b3 + b4 + b6 + b7 + b9 + b10 + b12)*(X^3 + X^4 + X^6 + X^7 + X^9 + X1 + X2 + X4 + X5 + X7 + X9 + X^20 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a3, a4, a6, a7, a9, a10, a12),
				mzd_sum(t2, 9, b0, b1, b3, b4, b6, b7, b9, b10, b12), 0), X, 14,  3, 4, 6, 7, 9, 11, 12, 14, 15, 17, 19, 20, 22, 23);
	/* (a0)(b0)*(1 + X + X^7 + X^8 + X5 + X6 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a0, b0, 0), X, 8,  0, 1, 7, 8, 15, 16, 22, 23);
	/* (a1)(b1)*(X + X^8 + X6 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a1, b1, 0), X, 4,  1, 8, 16, 23);
	/* (a0 + a1)(b0 + b1)*(X + X^8 + X6 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 4,  1, 8, 16, 23);
	/* (a0 + a2 + a4 + a6 + a8 + a10 + a12)(b0 + b2 + b4 + b6 + b8 + b10 + b12)*(X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X7 + X8 + X9 + X^20 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a4, a6, a8, a10, a12),
				mzd_sum(t2, 7, b0, b2, b4, b6, b8, b10, b12), 0), X, 14,  2, 3, 4, 5, 6, 7, 8, 17, 18, 19, 20, 21, 22, 23);
	/* (a1 + a3 + a5 + a7 + a9 + a11)(b1 + b3 + b5 + b7 + b9 + b11)*(X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X7 + X8 + X9 + X^20 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a5, a7, a9, a11),
				mzd_sum(t2, 6, b1, b3, b5, b7, b9, b11), 0), X, 14,  2, 3, 4, 5, 6, 7, 8, 17, 18, 19, 20, 21, 22, 23);
	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12)(b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10 + b11 + b12)*(X^3 + X^5 + X^7 + X^9 + X0 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X9 + X^21 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 13, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12),
				mzd_sum(t2, 13, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12), 0), X, 15,  3, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 21, 23);
	/* (a0 + a3 + a5 + a6 + a7 + a10 + a12)(b0 + b3 + b5 + b6 + b7 + b10 + b12)*(X^2 + X^3 + X^4 + X^6 + X7 + X8 + X9 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a3, a5, a6, a7, a10, a12),
				mzd_sum(t2, 7, b0, b3, b5, b6, b7, b10, b12), 0), X, 8,  2, 3, 4, 6, 17, 18, 19, 21);
	/* (a1 + a3 + a4 + a5 + a8 + a10 + a11 + a12)(b1 + b3 + b4 + b5 + b8 + b10 + b11 + b12)*(X^3 + X^4 + X^5 + X^7 + X8 + X9 + X^20 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a3, a4, a5, a8, a10, a11, a12),
				mzd_sum(t2, 8, b1, b3, b4, b5, b8, b10, b11, b12), 0), X, 8,  3, 4, 5, 7, 18, 19, 20, 22);
	/* (a2 + a4 + a5 + a6 + a9 + a11 + a12)(b2 + b4 + b5 + b6 + b9 + b11 + b12)*(X^4 + X^5 + X^6 + X^8 + X9 + X^20 + X^21 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a4, a5, a6, a9, a11, a12),
				mzd_sum(t2, 7, b2, b4, b5, b6, b9, b11, b12), 0), X, 8,  4, 5, 6, 8, 19, 20, 21, 23);
	/* (a0 + a1 + a4 + a6 + a7 + a8 + a11)(b0 + b1 + b4 + b6 + b7 + b8 + b11)*(X^2 + X^5 + X^6 + X^7 + X7 + X^20 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a4, a6, a7, a8, a11),
				mzd_sum(t2, 7, b0, b1, b4, b6, b7, b8, b11), 0), X, 8,  2, 5, 6, 7, 17, 20, 21, 22);
	/* (a0 + a2 + a3 + a4 + a7 + a9 + a10 + a11)(b0 + b2 + b3 + b4 + b7 + b9 + b10 + b11)*(X^3 + X^6 + X^7 + X^8 + X8 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a3, a4, a7, a9, a10, a11),
				mzd_sum(t2, 8, b0, b2, b3, b4, b7, b9, b10, b11), 0), X, 8,  3, 6, 7, 8, 18, 21, 22, 23);
	/* (a1 + a2 + a3 + a6 + a8 + a9 + a10)(b1 + b2 + b3 + b6 + b8 + b9 + b10)*(X^2 + X^4 + X^7 + X^8 + X7 + X9 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a2, a3, a6, a8, a9, a10),
				mzd_sum(t2, 7, b1, b2, b3, b6, b8, b9, b10), 0), X, 8,  2, 4, 7, 8, 17, 19, 22, 23);
	/* (a0 + a3 + a4 + a5 + a7 + a10 + a11 + a12)(b0 + b3 + b4 + b5 + b7 + b10 + b11 + b12)*(X^2 + X^4 + X^5 + X^6 + X7 + X9 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a3, a4, a5, a7, a10, a11, a12),
				mzd_sum(t2, 8, b0, b3, b4, b5, b7, b10, b11, b12), 0), X, 8,  2, 4, 5, 6, 17, 19, 20, 21);
	/* (a1 + a4 + a5 + a6 + a8 + a11 + a12)(b1 + b4 + b5 + b6 + b8 + b11 + b12)*(X^3 + X^5 + X^6 + X^7 + X8 + X^20 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a4, a5, a6, a8, a11, a12),
				mzd_sum(t2, 7, b1, b4, b5, b6, b8, b11, b12), 0), X, 8,  3, 5, 6, 7, 18, 20, 21, 22);
	/* (a2 + a3 + a4 + a6 + a9 + a10 + a11)(b2 + b3 + b4 + b6 + b9 + b10 + b11)*(X^4 + X^6 + X^7 + X^8 + X9 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a4, a6, a9, a10, a11),
				mzd_sum(t2, 7, b2, b3, b4, b6, b9, b10, b11), 0), X, 8,  4, 6, 7, 8, 19, 21, 22, 23);
	/* (a0 + a1 + a3 + a6 + a7 + a8 + a10)(b0 + b1 + b3 + b6 + b7 + b8 + b10)*(X^2 + X^3 + X^6 + X^8 + X7 + X8 + X^21 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a3, a6, a7, a8, a10),
				mzd_sum(t2, 7, b0, b1, b3, b6, b7, b8, b10), 0), X, 8,  2, 3, 6, 8, 17, 18, 21, 23);
	/* (a0 + a2 + a5 + a6 + a7 + a9 + a12)(b0 + b2 + b5 + b6 + b7 + b9 + b12)*(X^2 + X^3 + X^4 + X^7 + X7 + X8 + X9 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a5, a6, a7, a9, a12),
				mzd_sum(t2, 7, b0, b2, b5, b6, b7, b9, b12), 0), X, 8,  2, 3, 4, 7, 17, 18, 19, 22);
	/* (a1 + a2 + a3 + a5 + a8 + a9 + a10 + a12)(b1 + b2 + b3 + b5 + b8 + b9 + b10 + b12)*(X^3 + X^4 + X^5 + X^8 + X8 + X9 + X^20 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a2, a3, a5, a8, a9, a10, a12),
				mzd_sum(t2, 8, b1, b2, b3, b5, b8, b9, b10, b12), 0), X, 8,  3, 4, 5, 8, 18, 19, 20, 23);
	/* (a0 + a4 + a7 + a8 + a10 + a12)(b0 + b4 + b7 + b8 + b10 + b12)*(X^4 + X^5 + X^6 + X^7 + X^9 + X3 + X4 + X5 + X6 + X8 + X9 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a4, a7, a8, a10, a12),
				mzd_sum(t2, 6, b0, b4, b7, b8, b10, b12), 0), X, 12,  4, 5, 6, 7, 9, 13, 14, 15, 16, 18, 19, 22);
	/* (a1 + a4 + a5 + a7 + a9 + a10 + a11 + a12)(b1 + b4 + b5 + b7 + b9 + b10 + b11 + b12)*(X^5 + X^6 + X^7 + X^8 + X0 + X4 + X5 + X6 + X7 + X9 + X^20 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a4, a5, a7, a9, a10, a11, a12),
				mzd_sum(t2, 8, b1, b4, b5, b7, b9, b10, b11, b12), 0), X, 12,  5, 6, 7, 8, 10, 14, 15, 16, 17, 19, 20, 23);
	/* (a2 + a5 + a6 + a8 + a10 + a11 + a12)(b2 + b5 + b6 + b8 + b10 + b11 + b12)*(X^2 + X^6 + X^7 + X^8 + X1 + X5 + X6 + X8 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a5, a6, a8, a10, a11, a12),
				mzd_sum(t2, 7, b2, b5, b6, b8, b10, b11, b12), 0), X, 10,  2, 6, 7, 8, 11, 15, 16, 18, 20, 21);
	/* (a3 + a6 + a7 + a9 + a11 + a12)(b3 + b6 + b7 + b9 + b11 + b12)*(X^3 + X^7 + X^8 + X^9 + X2 + X6 + X7 + X9 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a3, a6, a7, a9, a11, a12),
				mzd_sum(t2, 6, b3, b6, b7, b9, b11, b12), 0), X, 10,  3, 7, 8, 9, 12, 16, 17, 19, 21, 22);
	/* (a0 + a1 + a2 + a3 + a7 + a10 + a11)(b0 + b1 + b2 + b3 + b7 + b10 + b11)*(X^2 + X^3 + X^6 + X^9 + X1 + X2 + X5 + X7 + X8 + X9 + X^20 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a2, a3, a7, a10, a11),
				mzd_sum(t2, 7, b0, b1, b2, b3, b7, b10, b11), 0), X, 12,  2, 3, 6, 9, 11, 12, 15, 17, 18, 19, 20, 22);
	/* (a0 + a1 + a5 + a8 + a9 + a11)(b0 + b1 + b5 + b8 + b9 + b11)*(X^2 + X^3 + X^4 + X^6 + X^8 + X0 + X1 + X2 + X3 + X5 + X9 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a5, a8, a9, a11),
				mzd_sum(t2, 6, b0, b1, b5, b8, b9, b11), 0), X, 12,  2, 3, 4, 6, 8, 10, 11, 12, 13, 15, 19, 23);
	/* (a0 + a2 + a4 + a5 + a6 + a7 + a11)(b0 + b2 + b4 + b5 + b6 + b7 + b11)*(X^3 + X^5 + X^6 + X^9 + X0 + X2 + X4 + X5 + X^20 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a4, a5, a6, a7, a11),
				mzd_sum(t2, 7, b0, b2, b4, b5, b6, b7, b11), 0), X, 12,  3, 5, 6, 9, 10, 12, 14, 15, 20, 21, 22, 23);
	/* (a1 + a3 + a4 + a5 + a6 + a10)(b1 + b3 + b4 + b5 + b6 + b10)*(X^2 + X^4 + X^6 + X^7 + X^9 + X0 + X1 + X3 + X5 + X6 + X7 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a4, a5, a6, a10),
				mzd_sum(t2, 6, b1, b3, b4, b5, b6, b10), 0), X, 14,  2, 4, 6, 7, 9, 10, 11, 13, 15, 16, 17, 21, 22, 23);
	/* (a2 + a3 + a5 + a7 + a8 + a9 + a10)(b2 + b3 + b5 + b7 + b8 + b9 + b10)*(X^3 + X^4 + X^5 + X^6 + X^8 + X2 + X3 + X4 + X5 + X7 + X8 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a5, a7, a8, a9, a10),
				mzd_sum(t2, 7, b2, b3, b5, b7, b8, b9, b10), 0), X, 12,  3, 4, 5, 6, 8, 12, 13, 14, 15, 17, 18, 21);
	/* (a0 + a4 + a5 + a6 + a7 + a9 + a11 + a12)(b0 + b4 + b5 + b6 + b7 + b9 + b11 + b12)*(X^2 + X^4 + X^5 + X^6 + X^7 + X^9 + X2 + X3 + X5 + X8 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a4, a5, a6, a7, a9, a11, a12),
				mzd_sum(t2, 8, b0, b4, b5, b6, b7, b9, b11, b12), 0), X, 12,  2, 4, 5, 6, 7, 9, 12, 13, 15, 18, 21, 22);
	/* (a1 + a5 + a6 + a7 + a8 + a10 + a12)(b1 + b5 + b6 + b7 + b8 + b10 + b12)*(X^3 + X^5 + X^6 + X^7 + X^8 + X0 + X3 + X4 + X6 + X9 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a5, a6, a7, a8, a10, a12),
				mzd_sum(t2, 7, b1, b5, b6, b7, b8, b10, b12), 0), X, 12,  3, 5, 6, 7, 8, 10, 13, 14, 16, 19, 22, 23);
	/* (a2 + a6 + a7 + a8 + a9 + a11)(b2 + b6 + b7 + b8 + b9 + b11)*(X^2 + X^4 + X^6 + X^7 + X^8 + X1 + X4 + X5 + X^20 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a6, a7, a8, a9, a11),
				mzd_sum(t2, 6, b2, b6, b7, b8, b9, b11), 0), X, 10,  2, 4, 6, 7, 8, 11, 14, 15, 20, 23);
	/* (a3 + a4 + a5 + a6 + a8 + a10 + a11)(b3 + b4 + b5 + b6 + b8 + b10 + b11)*(X^2 + X^3 + X^5 + X^7 + X^8 + X2 + X5 + X6 + X7 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a3, a4, a5, a6, a8, a10, a11),
				mzd_sum(t2, 7, b3, b4, b5, b6, b8, b10, b11), 0), X, 10,  2, 3, 5, 7, 8, 12, 15, 16, 17, 21);
	/* (a0 + a1 + a2 + a3 + a5 + a7 + a8 + a11)(b0 + b1 + b2 + b3 + b5 + b7 + b8 + b11)*(X^2 + X^3 + X^4 + X^8 + X^9 + X0 + X2 + X4 + X5 + X6 + X8 + X9 + X^21 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a2, a3, a5, a7, a8, a11),
				mzd_sum(t2, 8, b0, b1, b2, b3, b5, b7, b8, b11), 0), X, 14,  2, 3, 4, 8, 9, 10, 12, 14, 15, 16, 18, 19, 21, 23);
	/* (a0 + a1 + a4 + a8 + a9 + a10 + a11)(b0 + b1 + b4 + b8 + b9 + b10 + b11)*(X^3 + X^4 + X^6 + X^8 + X^9 + X3 + X6 + X7 + X8 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a4, a8, a9, a10, a11),
				mzd_sum(t2, 7, b0, b1, b4, b8, b9, b10, b11), 0), X, 10,  3, 4, 6, 8, 9, 13, 16, 17, 18, 22);
	/* (a0 + a2 + a4 + a5 + a8 + a12)(b0 + b2 + b4 + b5 + b8 + b12)*(X^4 + X^7 + X^8 + X0 + X1 + X2 + X3 + X7 + X9 + X^20 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a4, a5, a8, a12),
				mzd_sum(t2, 6, b0, b2, b4, b5, b8, b12), 0), X, 12,  4, 7, 8, 10, 11, 12, 13, 17, 19, 20, 21, 22);
	/* (a1 + a3 + a4 + a7 + a11 + a12)(b1 + b3 + b4 + b7 + b11 + b12)*(X^5 + X^8 + X^9 + X1 + X2 + X3 + X4 + X8 + X^20 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a4, a7, a11, a12),
				mzd_sum(t2, 6, b1, b3, b4, b7, b11, b12), 0), X, 12,  5, 8, 9, 11, 12, 13, 14, 18, 20, 21, 22, 23);
	/* (a2 + a3 + a4 + a5 + a7 + a9 + a10)(b2 + b3 + b4 + b5 + b7 + b9 + b10)*(X^2 + X^5 + X^6 + X^8 + X^9 + X0 + X1 + X5 + X7 + X8 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a4, a5, a7, a9, a10),
				mzd_sum(t2, 7, b2, b3, b4, b5, b7, b9, b10), 0), X, 12,  2, 5, 6, 8, 9, 10, 11, 15, 17, 18, 19, 20);
	/* (a0 + a4 + a5 + a9 + a10)(b0 + b4 + b5 + b9 + b10)*(X^5 + X^6 + X0 + X1 + X2 + X3 + X5 + X6 + X7 + X8 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a4, a5, a9, a10),
				mzd_sum(t2, 5, b0, b4, b5, b9, b10), 0), X, 12,  5, 6, 10, 11, 12, 13, 15, 16, 17, 18, 22, 23);
	/* (a1 + a4 + a6 + a9 + a11)(b1 + b4 + b6 + b9 + b11)*(X^2 + X^6 + X^7 + X^9 + X1 + X2 + X3 + X4 + X6 + X8 + X9 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a4, a6, a9, a11),
				mzd_sum(t2, 5, b1, b4, b6, b9, b11), 0), X, 12,  2, 6, 7, 9, 11, 12, 13, 14, 16, 18, 19, 23);
	/* (a2 + a4 + a7 + a9 + a12)(b2 + b4 + b7 + b9 + b12)*(X^2 + X^3 + X^7 + X^8 + X^9 + X0 + X2 + X3 + X4 + X5 + X9 + X^20) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a4, a7, a9, a12),
				mzd_sum(t2, 5, b2, b4, b7, b9, b12), 0), X, 12,  2, 3, 7, 8, 9, 10, 12, 13, 14, 15, 19, 20);
	/* (a3 + a4 + a8 + a9)(b3 + b4 + b8 + b9)*(X^3 + X^4 + X^8 + X^9 + X0 + X1 + X3 + X4 + X5 + X6 + X^20 + X^21) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a3, a4, a8, a9),
				mzd_sum(t2, 4, b3, b4, b8, b9), 0), X, 12,  3, 4, 8, 9, 10, 11, 13, 14, 15, 16, 20, 21);
	/* (a0 + a1 + a2 + a3 + a5 + a6 + a7 + a8 + a10 + a11 + a12)(b0 + b1 + b2 + b3 + b5 + b6 + b7 + b8 + b10 + b11 + b12)*(X^4 + X^5 + X^9 + X0 + X1 + X2 + X4 + X5 + X6 + X7 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 11, a0, a1, a2, a3, a5, a6, a7, a8, a10, a11, a12),
				mzd_sum(t2, 11, b0, b1, b2, b3, b5, b6, b7, b8, b10, b11, b12), 0), X, 12,  4, 5, 9, 10, 11, 12, 14, 15, 16, 17, 21, 22);
	/* (a0 + a1 + a5 + a6 + a10 + a11)(b0 + b1 + b5 + b6 + b10 + b11)*(X^2 + X^3 + X^4 + X^5 + X^7 + X^8 + X1 + X3 + X6 + X7 + X9 + X^20 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a5, a6, a10, a11),
				mzd_sum(t2, 6, b0, b1, b5, b6, b10, b11), 0), X, 14,  2, 3, 4, 5, 7, 8, 11, 13, 16, 17, 19, 20, 21, 22);
	/* (a0 + a2 + a5 + a7 + a10 + a12)(b0 + b2 + b5 + b7 + b10 + b12)*(X^3 + X^5 + X^8 + X2 + X3 + X7 + X^20 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a5, a7, a10, a12),
				mzd_sum(t2, 6, b0, b2, b5, b7, b10, b12), 0), X, 8,  3, 5, 8, 12, 13, 17, 20, 22);
	/* (a1 + a3 + a6 + a8 + a11)(b1 + b3 + b6 + b8 + b11)*(X^4 + X^6 + X^9 + X3 + X4 + X8 + X^21 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a3, a6, a8, a11),
				mzd_sum(t2, 5, b1, b3, b6, b8, b11), 0), X, 8,  4, 6, 9, 13, 14, 18, 21, 23);
	/* (a2 + a3 + a7 + a8 + a12)(b2 + b3 + b7 + b8 + b12)*(X^2 + X^4 + X^5 + X^6 + X^7 + X0 + X3 + X5 + X7 + X8 + X9 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a3, a7, a8, a12),
				mzd_sum(t2, 5, b2, b3, b7, b8, b12), 0), X, 14,  2, 4, 5, 6, 7, 10, 13, 15, 17, 18, 19, 21, 22, 23);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}
void _mzd_ptr_addmul_karatsuba14(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 55 multiplications and three temporaries (not optimal!) */
	const mzd_t *a0  = A[ 0];
	const mzd_t *a1  = A[ 1];
	const mzd_t *a2  = A[ 2];
	const mzd_t *a3  = A[ 3];
	const mzd_t *a4  = A[ 4];
	const mzd_t *a5  = A[ 5];
	const mzd_t *a6  = A[ 6];
	const mzd_t *a7  = A[ 7];
	const mzd_t *a8  = A[ 8];
	const mzd_t *a9  = A[ 9];
	const mzd_t *a10 = A[10];
	const mzd_t *a11 = A[11];
	const mzd_t *a12 = A[12];
	const mzd_t *a13 = A[13];

	const mzd_t *b0  = B[ 0];
	const mzd_t *b1  = B[ 1];
	const mzd_t *b2  = B[ 2];
	const mzd_t *b3  = B[ 3];
	const mzd_t *b4  = B[ 4];
	const mzd_t *b5  = B[ 5];
	const mzd_t *b6  = B[ 6];
	const mzd_t *b7  = B[ 7];
	const mzd_t *b8  = B[ 8];
	const mzd_t *b9  = B[ 9];
	const mzd_t *b10 = B[10];
	const mzd_t *b11 = B[11];
	const mzd_t *b12 = B[12];
	const mzd_t *b13 = B[13];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	/* (a13)(b13)*(X^2 + X^4 + X^9 + X1 + X7 + X9 + X^24 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a13, b13, 0), X, 8,  2, 4, 9, 11, 17, 19, 24, 26);
	/* (a0 + a2 + a3 + a5 + a6 + a8 + a9 + a11 + a12)(b0 + b2 + b3 + b5 + b6 + b8 + b9 + b11 + b12)*(X^3 + X^4 + X^5 + X^7 + X^8 + X2 + X3 + X5 + X6 + X^20 + X^21 + X^23 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a3, a5, a6, a8, a9, a11, a12),
				mzd_sum(t2, 9, b0, b2, b3, b5, b6, b8, b9, b11, b12), 0), X, 14,  3, 4, 5, 7, 8, 12, 13, 15, 16, 20, 21, 23, 24, 25);
	/* (a1 + a2 + a4 + a5 + a7 + a8 + a10 + a11 + a13)(b1 + b2 + b4 + b5 + b7 + b8 + b10 + b11 + b13)*(X^2 + X^5 + X^6 + X^8 + X1 + X3 + X4 + X6 + X9 + X^21 + X^22 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a1, a2, a4, a5, a7, a8, a10, a11, a13),
				mzd_sum(t2, 9, b1, b2, b4, b5, b7, b8, b10, b11, b13), 0), X, 12,  2, 5, 6, 8, 11, 13, 14, 16, 19, 21, 22, 25);
	/* (a0 + a1 + a3 + a4 + a6 + a7 + a9 + a10 + a12 + a13)(b0 + b1 + b3 + b4 + b6 + b7 + b9 + b10 + b12 + b13)*(X^2 + X^3 + X^4 + X^6 + X^7 + X1 + X2 + X4 + X5 + X9 + X^20 + X^22 + X^23 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a1, a3, a4, a6, a7, a9, a10, a12, a13),
				mzd_sum(t2, 10, b0, b1, b3, b4, b6, b7, b9, b10, b12, b13), 0), X, 14,  2, 3, 4, 6, 7, 11, 12, 14, 15, 19, 20, 22, 23, 24);
	/* (a0)(b0)*(1 + X + X^2 + X^3 + X^7 + X^8 + X^9 + X0 + X5 + X6 + X7 + X8 + X^22 + X^23 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a0, b0, 0), X, 16,  0, 1, 2, 3, 7, 8, 9, 10, 15, 16, 17, 18, 22, 23, 24, 25);
	/* (a1)(b1)*(X + X^3 + X^8 + X0 + X6 + X8 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a1, b1, 0), X, 8,  1, 3, 8, 10, 16, 18, 23, 25);
	/* (a0 + a1)(b0 + b1)*(X + X^3 + X^8 + X0 + X6 + X8 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 8,  1, 3, 8, 10, 16, 18, 23, 25);
	/* (a0 + a3 + a5 + a6 + a7 + a10 + a12 + a13)(b0 + b3 + b5 + b6 + b7 + b10 + b12 + b13)*(X^3 + X^4 + X^6 + X^9 + X8 + X9 + X^21 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a3, a5, a6, a7, a10, a12, a13),
				mzd_sum(t2, 8, b0, b3, b5, b6, b7, b10, b12, b13), 0), X, 8,  3, 4, 6, 9, 18, 19, 21, 24);
	/* (a1 + a3 + a4 + a5 + a8 + a10 + a11 + a12)(b1 + b3 + b4 + b5 + b8 + b10 + b11 + b12)*(X^4 + X^5 + X^7 + X0 + X9 + X^20 + X^22 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a3, a4, a5, a8, a10, a11, a12),
				mzd_sum(t2, 8, b1, b3, b4, b5, b8, b10, b11, b12), 0), X, 8,  4, 5, 7, 10, 19, 20, 22, 25);
	/* (a2 + a4 + a5 + a6 + a9 + a11 + a12 + a13)(b2 + b4 + b5 + b6 + b9 + b11 + b12 + b13)*(X^2 + X^4 + X^5 + X^6 + X^8 + X^9 + X7 + X9 + X^20 + X^21 + X^23 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a4, a5, a6, a9, a11, a12, a13),
				mzd_sum(t2, 8, b2, b4, b5, b6, b9, b11, b12, b13), 0), X, 12,  2, 4, 5, 6, 8, 9, 17, 19, 20, 21, 23, 24);
	/* (a0 + a1 + a4 + a6 + a7 + a8 + a11 + a13)(b0 + b1 + b4 + b6 + b7 + b8 + b11 + b13)*(X^3 + X^5 + X^6 + X^7 + X^9 + X0 + X8 + X^20 + X^21 + X^22 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a4, a6, a7, a8, a11, a13),
				mzd_sum(t2, 8, b0, b1, b4, b6, b7, b8, b11, b13), 0), X, 12,  3, 5, 6, 7, 9, 10, 18, 20, 21, 22, 24, 25);
	/* (a0 + a2 + a3 + a4 + a7 + a9 + a10 + a11)(b0 + b2 + b3 + b4 + b7 + b9 + b10 + b11)*(X^2 + X^6 + X^7 + X^8 + X^9 + X0 + X7 + X^21 + X^22 + X^23 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a3, a4, a7, a9, a10, a11),
				mzd_sum(t2, 8, b0, b2, b3, b4, b7, b9, b10, b11), 0), X, 12,  2, 6, 7, 8, 9, 10, 17, 21, 22, 23, 24, 25);
	/* (a1 + a2 + a3 + a6 + a8 + a9 + a10 + a13)(b1 + b2 + b3 + b6 + b8 + b9 + b10 + b13)*(X^2 + X^3 + X^4 + X^7 + X^8 + X0 + X7 + X8 + X9 + X^22 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a2, a3, a6, a8, a9, a10, a13),
				mzd_sum(t2, 8, b1, b2, b3, b6, b8, b9, b10, b13), 0), X, 12,  2, 3, 4, 7, 8, 10, 17, 18, 19, 22, 23, 25);
	/* (a0 + a3 + a4 + a5 + a7 + a10 + a11 + a12)(b0 + b3 + b4 + b5 + b7 + b10 + b11 + b12)*(X^2 + X^3 + X^4 + X^5 + X^6 + X0 + X7 + X8 + X9 + X^20 + X^21 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a3, a4, a5, a7, a10, a11, a12),
				mzd_sum(t2, 8, b0, b3, b4, b5, b7, b10, b11, b12), 0), X, 12,  2, 3, 4, 5, 6, 10, 17, 18, 19, 20, 21, 25);
	/* (a1 + a4 + a5 + a6 + a8 + a11 + a12 + a13)(b1 + b4 + b5 + b6 + b8 + b11 + b12 + b13)*(X^2 + X^3 + X^5 + X^6 + X^7 + X^9 + X7 + X8 + X^20 + X^21 + X^22 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a4, a5, a6, a8, a11, a12, a13),
				mzd_sum(t2, 8, b1, b4, b5, b6, b8, b11, b12, b13), 0), X, 12,  2, 3, 5, 6, 7, 9, 17, 18, 20, 21, 22, 24);
	/* (a2 + a3 + a4 + a6 + a9 + a10 + a11 + a13)(b2 + b3 + b4 + b6 + b9 + b10 + b11 + b13)*(X^3 + X^4 + X^6 + X^7 + X^8 + X0 + X8 + X9 + X^21 + X^22 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a3, a4, a6, a9, a10, a11, a13),
				mzd_sum(t2, 8, b2, b3, b4, b6, b9, b10, b11, b13), 0), X, 12,  3, 4, 6, 7, 8, 10, 18, 19, 21, 22, 23, 25);
	/* (a0 + a1 + a3 + a6 + a7 + a8 + a10 + a13)(b0 + b1 + b3 + b6 + b7 + b8 + b10 + b13)*(X^3 + X^6 + X^8 + X^9 + X8 + X^21 + X^23 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a3, a6, a7, a8, a10, a13),
				mzd_sum(t2, 8, b0, b1, b3, b6, b7, b8, b10, b13), 0), X, 8,  3, 6, 8, 9, 18, 21, 23, 24);
	/* (a0 + a2 + a5 + a6 + a7 + a9 + a12 + a13)(b0 + b2 + b5 + b6 + b7 + b9 + b12 + b13)*(X^4 + X^7 + X^9 + X0 + X9 + X^22 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a5, a6, a7, a9, a12, a13),
				mzd_sum(t2, 8, b0, b2, b5, b6, b7, b9, b12, b13), 0), X, 8,  4, 7, 9, 10, 19, 22, 24, 25);
	/* (a1 + a2 + a3 + a5 + a8 + a9 + a10 + a12)(b1 + b2 + b3 + b5 + b8 + b9 + b10 + b12)*(X^2 + X^4 + X^5 + X^8 + X^9 + X0 + X7 + X9 + X^20 + X^23 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a2, a3, a5, a8, a9, a10, a12),
				mzd_sum(t2, 8, b1, b2, b3, b5, b8, b9, b10, b12), 0), X, 12,  2, 4, 5, 8, 9, 10, 17, 19, 20, 23, 24, 25);
	/* (a0 + a4 + a7 + a8 + a10 + a12 + a13)(b0 + b4 + b7 + b8 + b10 + b12 + b13)*(X^3 + X^4 + X^5 + X^6 + X^7 + X^9 + X0 + X3 + X4 + X5 + X6 + X9 + X^22 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a4, a7, a8, a10, a12, a13),
				mzd_sum(t2, 7, b0, b4, b7, b8, b10, b12, b13), 0), X, 14,  3, 4, 5, 6, 7, 9, 10, 13, 14, 15, 16, 19, 22, 25);
	/* (a1 + a4 + a5 + a7 + a9 + a10 + a11 + a12)(b1 + b4 + b5 + b7 + b9 + b10 + b11 + b12)*(X^2 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X4 + X5 + X6 + X9 + X^20 + X^23 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a4, a5, a7, a9, a10, a11, a12),
				mzd_sum(t2, 8, b1, b4, b5, b7, b9, b10, b11, b12), 0), X, 14,  2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 19, 20, 23, 24);
	/* (a2 + a5 + a6 + a8 + a10 + a11 + a12 + a13)(b2 + b5 + b6 + b8 + b10 + b11 + b12 + b13)*(X^3 + X^6 + X^7 + X^8 + X^9 + X0 + X1 + X5 + X6 + X7 + X^20 + X^21 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a5, a6, a8, a10, a11, a12, a13),
				mzd_sum(t2, 8, b2, b5, b6, b8, b10, b11, b12, b13), 0), X, 14,  3, 6, 7, 8, 9, 10, 11, 15, 16, 17, 20, 21, 24, 25);
	/* (a3 + a6 + a7 + a9 + a11 + a12 + a13)(b3 + b6 + b7 + b9 + b11 + b12 + b13)*(X^2 + X^7 + X^8 + X0 + X2 + X6 + X8 + X9 + X^21 + X^22 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a3, a6, a7, a9, a11, a12, a13),
				mzd_sum(t2, 7, b3, b6, b7, b9, b11, b12, b13), 0), X, 12,  2, 7, 8, 10, 12, 16, 18, 19, 21, 22, 24, 25);
	/* (a0 + a1 + a2 + a3 + a7 + a10 + a11 + a13)(b0 + b1 + b2 + b3 + b7 + b10 + b11 + b13)*(X^2 + X^3 + X^6 + X^9 + X1 + X2 + X5 + X7 + X8 + X9 + X^20 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a2, a3, a7, a10, a11, a13),
				mzd_sum(t2, 8, b0, b1, b2, b3, b7, b10, b11, b13), 0), X, 12,  2, 3, 6, 9, 11, 12, 15, 17, 18, 19, 20, 22);
	/* (a0 + a1 + a5 + a8 + a9 + a11 + a13)(b0 + b1 + b5 + b8 + b9 + b11 + b13)*(X^4 + X^6 + X^8 + X^9 + X1 + X2 + X3 + X5 + X7 + X8 + X9 + X^23 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a5, a8, a9, a11, a13),
				mzd_sum(t2, 7, b0, b1, b5, b8, b9, b11, b13), 0), X, 14,  4, 6, 8, 9, 11, 12, 13, 15, 17, 18, 19, 23, 24, 25);
	/* (a0 + a2 + a4 + a5 + a6 + a7 + a11)(b0 + b2 + b4 + b5 + b6 + b7 + b11)*(X^5 + X^6 + X^9 + X2 + X4 + X5 + X8 + X^20 + X^21 + X^22 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a4, a5, a6, a7, a11),
				mzd_sum(t2, 7, b0, b2, b4, b5, b6, b7, b11), 0), X, 12,  5, 6, 9, 12, 14, 15, 18, 20, 21, 22, 23, 25);
	/* (a1 + a3 + a4 + a5 + a6 + a10 + a13)(b1 + b3 + b4 + b5 + b6 + b10 + b13)*(X^2 + X^4 + X^6 + X^7 + X^9 + X0 + X1 + X3 + X5 + X6 + X7 + X^21 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a3, a4, a5, a6, a10, a13),
				mzd_sum(t2, 7, b1, b3, b4, b5, b6, b10, b13), 0), X, 14,  2, 4, 6, 7, 9, 10, 11, 13, 15, 16, 17, 21, 22, 23);
	/* (a2 + a3 + a5 + a7 + a8 + a9 + a10)(b2 + b3 + b5 + b7 + b8 + b9 + b10)*(X^2 + X^3 + X^4 + X^5 + X^6 + X^8 + X^9 + X2 + X3 + X4 + X5 + X8 + X^21 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a3, a5, a7, a8, a9, a10),
				mzd_sum(t2, 7, b2, b3, b5, b7, b8, b9, b10), 0), X, 14,  2, 3, 4, 5, 6, 8, 9, 12, 13, 14, 15, 18, 21, 24);
	/* (a0 + a4 + a5 + a6 + a7 + a9 + a11 + a12)(b0 + b4 + b5 + b6 + b7 + b9 + b11 + b12)*(X^2 + X^4 + X^5 + X^6 + X^7 + X^9 + X2 + X3 + X5 + X8 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a4, a5, a6, a7, a9, a11, a12),
				mzd_sum(t2, 8, b0, b4, b5, b6, b7, b9, b11, b12), 0), X, 12,  2, 4, 5, 6, 7, 9, 12, 13, 15, 18, 21, 22);
	/* (a1 + a5 + a6 + a7 + a8 + a10 + a12 + a13)(b1 + b5 + b6 + b7 + b8 + b10 + b12 + b13)*(X^3 + X^5 + X^6 + X^7 + X^8 + X0 + X3 + X4 + X6 + X9 + X^22 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a5, a6, a7, a8, a10, a12, a13),
				mzd_sum(t2, 8, b1, b5, b6, b7, b8, b10, b12, b13), 0), X, 12,  3, 5, 6, 7, 8, 10, 13, 14, 16, 19, 22, 23);
	/* (a2 + a6 + a7 + a8 + a9 + a11 + a13)(b2 + b6 + b7 + b8 + b9 + b11 + b13)*(X^4 + X^6 + X^7 + X^8 + X^9 + X1 + X4 + X5 + X7 + X^20 + X^23 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a2, a6, a7, a8, a9, a11, a13),
				mzd_sum(t2, 7, b2, b6, b7, b8, b9, b11, b13), 0), X, 12,  4, 6, 7, 8, 9, 11, 14, 15, 17, 20, 23, 24);
	/* (a3 + a4 + a5 + a6 + a8 + a10 + a11)(b3 + b4 + b5 + b6 + b8 + b10 + b11)*(X^5 + X^7 + X^8 + X^9 + X0 + X2 + X5 + X6 + X8 + X^21 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a3, a4, a5, a6, a8, a10, a11),
				mzd_sum(t2, 7, b3, b4, b5, b6, b8, b10, b11), 0), X, 12,  5, 7, 8, 9, 10, 12, 15, 16, 18, 21, 24, 25);
	/* (a0 + a1 + a2 + a3 + a5 + a7 + a8 + a11)(b0 + b1 + b2 + b3 + b5 + b7 + b8 + b11)*(X^2 + X^3 + X^4 + X^8 + X^9 + X0 + X2 + X4 + X5 + X6 + X8 + X9 + X^21 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a2, a3, a5, a7, a8, a11),
				mzd_sum(t2, 8, b0, b1, b2, b3, b5, b7, b8, b11), 0), X, 14,  2, 3, 4, 8, 9, 10, 12, 14, 15, 16, 18, 19, 21, 23);
	/* (a0 + a1 + a4 + a8 + a9 + a10 + a11 + a13)(b0 + b1 + b4 + b8 + b9 + b10 + b11 + b13)*(X^2 + X^4 + X^6 + X^8 + X0 + X3 + X6 + X^22 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a4, a8, a9, a10, a11, a13),
				mzd_sum(t2, 8, b0, b1, b4, b8, b9, b10, b11, b13), 0), X, 10,  2, 4, 6, 8, 10, 13, 16, 22, 24, 25);
	/* (a0 + a2 + a4 + a5 + a8 + a12 + a13)(b0 + b2 + b4 + b5 + b8 + b12 + b13)*(X^3 + X^4 + X^7 + X^8 + X1 + X2 + X3 + X7 + X8 + X9 + X^20 + X^21 + X^22 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a4, a5, a8, a12, a13),
				mzd_sum(t2, 7, b0, b2, b4, b5, b8, b12, b13), 0), X, 14,  3, 4, 7, 8, 11, 12, 13, 17, 18, 19, 20, 21, 22, 25);
	/* (a1 + a3 + a4 + a7 + a11 + a12 + a13)(b1 + b3 + b4 + b7 + b11 + b12 + b13)*(X^2 + X^5 + X^8 + X1 + X2 + X3 + X4 + X7 + X8 + X^20 + X^21 + X^22 + X^23 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a3, a4, a7, a11, a12, a13),
				mzd_sum(t2, 7, b1, b3, b4, b7, b11, b12, b13), 0), X, 14,  2, 5, 8, 11, 12, 13, 14, 17, 18, 20, 21, 22, 23, 24);
	/* (a2 + a3 + a4 + a5 + a7 + a9 + a10 + a13)(b2 + b3 + b4 + b5 + b7 + b9 + b10 + b13)*(X^2 + X^3 + X^5 + X^6 + X^8 + X^9 + X1 + X5 + X7 + X9 + X^20 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a3, a4, a5, a7, a9, a10, a13),
				mzd_sum(t2, 8, b2, b3, b4, b5, b7, b9, b10, b13), 0), X, 12,  2, 3, 5, 6, 8, 9, 11, 15, 17, 19, 20, 25);
	/* (a0 + a4 + a5 + a9 + a10)(b0 + b4 + b5 + b9 + b10)*(X^3 + X^5 + X^6 + X1 + X2 + X3 + X5 + X6 + X7 + X^22 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a0, a4, a5, a9, a10),
				mzd_sum(t2, 5, b0, b4, b5, b9, b10), 0), X, 12,  3, 5, 6, 11, 12, 13, 15, 16, 17, 22, 23, 25);
	/* (a1 + a4 + a6 + a9 + a11)(b1 + b4 + b6 + b9 + b11)*(X^2 + X^6 + X^7 + X^9 + X1 + X2 + X3 + X4 + X6 + X8 + X9 + X^23) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a4, a6, a9, a11),
				mzd_sum(t2, 5, b1, b4, b6, b9, b11), 0), X, 12,  2, 6, 7, 9, 11, 12, 13, 14, 16, 18, 19, 23);
	/* (a2 + a4 + a7 + a9 + a12)(b2 + b4 + b7 + b9 + b12)*(X^3 + X^7 + X^8 + X0 + X2 + X3 + X4 + X5 + X7 + X9 + X^20 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a2, a4, a7, a9, a12),
				mzd_sum(t2, 5, b2, b4, b7, b9, b12), 0), X, 12,  3, 7, 8, 10, 12, 13, 14, 15, 17, 19, 20, 24);
	/* (a3 + a4 + a8 + a9 + a13)(b3 + b4 + b8 + b9 + b13)*(X^4 + X^8 + X^9 + X1 + X3 + X4 + X5 + X6 + X8 + X^20 + X^21 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a3, a4, a8, a9, a13),
				mzd_sum(t2, 5, b3, b4, b8, b9, b13), 0), X, 12,  4, 8, 9, 11, 13, 14, 15, 16, 18, 20, 21, 25);
	/* (a0 + a1 + a2 + a3 + a5 + a6 + a7 + a8 + a10 + a11 + a12 + a13)(b0 + b1 + b2 + b3 + b5 + b6 + b7 + b8 + b10 + b11 + b12 + b13)*(X^2 + X^4 + X^5 + X0 + X1 + X2 + X4 + X5 + X6 + X^21 + X^22 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 12, a0, a1, a2, a3, a5, a6, a7, a8, a10, a11, a12, a13),
				mzd_sum(t2, 12, b0, b1, b2, b3, b5, b6, b7, b8, b10, b11, b12, b13), 0), X, 12,  2, 4, 5, 10, 11, 12, 14, 15, 16, 21, 22, 24);
	/* (a0 + a1 + a5 + a6 + a10 + a11)(b0 + b1 + b5 + b6 + b10 + b11)*(X^2 + X^3 + X^4 + X^5 + X^7 + X^8 + X1 + X3 + X6 + X7 + X9 + X^20 + X^21 + X^22) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a5, a6, a10, a11),
				mzd_sum(t2, 6, b0, b1, b5, b6, b10, b11), 0), X, 14,  2, 3, 4, 5, 7, 8, 11, 13, 16, 17, 19, 20, 21, 22);
	/* (a0 + a2 + a5 + a7 + a10 + a12)(b0 + b2 + b5 + b7 + b10 + b12)*(X^2 + X^5 + X^8 + X^9 + X0 + X2 + X3 + X8 + X^20 + X^22 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a5, a7, a10, a12),
				mzd_sum(t2, 6, b0, b2, b5, b7, b10, b12), 0), X, 12,  2, 5, 8, 9, 10, 12, 13, 18, 20, 22, 24, 25);
	/* (a1 + a3 + a6 + a8 + a11 + a13)(b1 + b3 + b6 + b8 + b11 + b13)*(X^2 + X^3 + X^4 + X^6 + X0 + X3 + X4 + X7 + X^21 + X^23 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a6, a8, a11, a13),
				mzd_sum(t2, 6, b1, b3, b6, b8, b11, b13), 0), X, 12,  2, 3, 4, 6, 10, 13, 14, 17, 21, 23, 24, 25);
	/* (a2 + a3 + a7 + a8 + a12 + a13)(b2 + b3 + b7 + b8 + b12 + b13)*(X^4 + X^5 + X^6 + X^7 + X^9 + X0 + X3 + X5 + X8 + X9 + X^21 + X^22 + X^23 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a3, a7, a8, a12, a13),
				mzd_sum(t2, 6, b2, b3, b7, b8, b12, b13), 0), X, 14,  4, 5, 6, 7, 9, 10, 13, 15, 18, 19, 21, 22, 23, 24);
	/* (a0 + a4 + a8 + a12)(b0 + b4 + b8 + b12)*(X^2 + X^3 + X^9 + X0 + X7 + X8 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a4, a8, a12),
				mzd_sum(t2, 4, b0, b4, b8, b12), 0), X, 8,  2, 3, 9, 10, 17, 18, 24, 25);
	/* (a1 + a5 + a9 + a13)(b1 + b5 + b9 + b13)*(X^2 + X^3 + X^9 + X0 + X7 + X8 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a1, a5, a9, a13),
				mzd_sum(t2, 4, b1, b5, b9, b13), 0), X, 8,  2, 3, 9, 10, 17, 18, 24, 25);
	/* (a2 + a6 + a10)(b2 + b6 + b10)*(X^2 + X^3 + X^9 + X0 + X7 + X8 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 3, a2, a6, a10),
				mzd_sum(t2, 3, b2, b6, b10), 0), X, 8,  2, 3, 9, 10, 17, 18, 24, 25);
	/* (a3 + a7 + a11)(b3 + b7 + b11)*(X^2 + X^3 + X^9 + X0 + X7 + X8 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 3, a3, a7, a11),
				mzd_sum(t2, 3, b3, b7, b11), 0), X, 8,  2, 3, 9, 10, 17, 18, 24, 25);
	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13)(b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10 + b11 + b12 + b13)*(X^5 + X^7 + X^9 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X^21 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 14, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13),
				mzd_sum(t2, 14, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13), 0), X, 15,  5, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 23, 25);
	/* (a0 + a1 + a4 + a5 + a8 + a9 + a12 + a13)(b0 + b1 + b4 + b5 + b8 + b9 + b12 + b13)*(X^3 + X0 + X8 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a4, a5, a8, a9, a12, a13),
				mzd_sum(t2, 8, b0, b1, b4, b5, b8, b9, b12, b13), 0), X, 4,  3, 10, 18, 25);
	/* (a0 + a2 + a4 + a6 + a8 + a10 + a12)(b0 + b2 + b4 + b6 + b8 + b10 + b12)*(X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X9 + X^20 + X^21 + X^22 + X^23 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a4, a6, a8, a10, a12),
				mzd_sum(t2, 7, b0, b2, b4, b6, b8, b10, b12), 0), X, 14,  4, 5, 6, 7, 8, 9, 10, 19, 20, 21, 22, 23, 24, 25);
	/* (a1 + a3 + a5 + a7 + a9 + a11 + a13)(b1 + b3 + b5 + b7 + b9 + b11 + b13)*(X^2 + X^4 + X^5 + X^6 + X^7 + X^8 + X0 + X7 + X9 + X^20 + X^21 + X^22 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a3, a5, a7, a9, a11, a13),
				mzd_sum(t2, 7, b1, b3, b5, b7, b9, b11, b13), 0), X, 14,  2, 4, 5, 6, 7, 8, 10, 17, 19, 20, 21, 22, 23, 25);
	/* (a2 + a3 + a6 + a7 + a10 + a11)(b2 + b3 + b6 + b7 + b10 + b11)*(X^3 + X0 + X8 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a3, a6, a7, a10, a11),
				mzd_sum(t2, 6, b2, b3, b6, b7, b10, b11), 0), X, 4,  3, 10, 18, 25);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba15(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 60 multiplications and three temporaries (not optimal!) */
	const mzd_t *a0  = A[ 0];
	const mzd_t *a1  = A[ 1];
	const mzd_t *a2  = A[ 2];
	const mzd_t *a3  = A[ 3];
	const mzd_t *a4  = A[ 4];
	const mzd_t *a5  = A[ 5];
	const mzd_t *a6  = A[ 6];
	const mzd_t *a7  = A[ 7];
	const mzd_t *a8  = A[ 8];
	const mzd_t *a9  = A[ 9];
	const mzd_t *a10 = A[10];
	const mzd_t *a11 = A[11];
	const mzd_t *a12 = A[12];
	const mzd_t *a13 = A[13];
	const mzd_t *a14 = A[14];

	const mzd_t *b0  = B[ 0];
	const mzd_t *b1  = B[ 1];
	const mzd_t *b2  = B[ 2];
	const mzd_t *b3  = B[ 3];
	const mzd_t *b4  = B[ 4];
	const mzd_t *b5  = B[ 5];
	const mzd_t *b6  = B[ 6];
	const mzd_t *b7  = B[ 7];
	const mzd_t *b8  = B[ 8];
	const mzd_t *b9  = B[ 9];
	const mzd_t *b10 = B[10];
	const mzd_t *b11 = B[11];
	const mzd_t *b12 = B[12];
	const mzd_t *b13 = B[13];
	const mzd_t *b14 = B[14];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	/* (a14)(b14)*(X + X^3 + X^6 + X^8 + X0 + X3 + X6 + X8 + X^21 + X^23 + X^25 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a14, b14, 0), X, 12,  1, 3, 6, 8, 10, 13, 16, 18, 21, 23, 25, 28);
	/* (a0)(b0)*(1 + X^2 + X^5 + X^7 + X^9 + X2 + X5 + X7 + X^20 + X^22 + X^24 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a0, b0, 0), X, 12,  0, 2, 5, 7, 9, 12, 15, 17, 20, 22, 24, 27);
	/* (a0 + a2 + a3 + a5 + a6 + a8 + a9 + a11 + a12 + a14)(b0 + b2 + b3 + b5 + b6 + b8 + b9 + b11 + b12 + b14)*(X + X^2 + X^3 + X^5 + X^7 + X^9 + X1 + X2 + X3 + X5 + X7 + X9 + X^20 + X^21 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a2, a3, a5, a6, a8, a9, a11, a12, a14),
				mzd_sum(t2, 10, b0, b2, b3, b5, b6, b8, b9, b11, b12, b14), 0), X, 16,  1, 2, 3, 5, 7, 9, 11, 12, 13, 15, 17, 19, 20, 21, 25, 26);
	/* (a1 + a2 + a4 + a5 + a7 + a8 + a10 + a11 + a13 + a14)(b1 + b2 + b4 + b5 + b7 + b8 + b10 + b11 + b13 + b14)*(X^2 + X^3 + X^4 + X^6 + X^8 + X0 + X2 + X3 + X4 + X6 + X8 + X^20 + X^21 + X^22 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a1, a2, a4, a5, a7, a8, a10, a11, a13, a14),
				mzd_sum(t2, 10, b1, b2, b4, b5, b7, b8, b10, b11, b13, b14), 0), X, 16,  2, 3, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 21, 22, 26, 27);
	/* (a0 + a1 + a3 + a4 + a6 + a7 + a9 + a10 + a12 + a13)(b0 + b1 + b3 + b4 + b6 + b7 + b9 + b10 + b12 + b13)*(X + X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X1 + X4 + X5 + X6 + X7 + X8 + X9 + X^22 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a1, a3, a4, a6, a7, a9, a10, a12, a13),
				mzd_sum(t2, 10, b0, b1, b3, b4, b6, b7, b9, b10, b12, b13), 0), X, 18,  1, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 22, 25, 27);
	/* (a0 + a2 + a4 + a6 + a8 + a10 + a12 + a14)(b0 + b2 + b4 + b6 + b8 + b10 + b12 + b14)*(X + X^2 + X^6 + X^7 + X0 + X1 + X2 + X6 + X7 + X^21 + X^22 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a4, a6, a8, a10, a12, a14),
				mzd_sum(t2, 8, b0, b2, b4, b6, b8, b10, b12, b14), 0), X, 14,  1, 2, 6, 7, 10, 11, 12, 16, 17, 21, 22, 25, 26, 27);
	/* (a1 + a3 + a5 + a7 + a9 + a11 + a13)(b1 + b3 + b5 + b7 + b9 + b11 + b13)*(X + X^2 + X^6 + X^7 + X0 + X1 + X2 + X6 + X7 + X^21 + X^22 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a3, a5, a7, a9, a11, a13),
				mzd_sum(t2, 7, b1, b3, b5, b7, b9, b11, b13), 0), X, 14,  1, 2, 6, 7, 10, 11, 12, 16, 17, 21, 22, 25, 26, 27);
	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14)(b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10 + b11 + b12 + b13 + b14)*(X^2 + X^7 + X1 + X3 + X4 + X5 + X6 + X8 + X9 + X^20 + X^21 + X^23 + X^24 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 15, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14),
				mzd_sum(t2, 15, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14), 0), X, 15,  2, 7, 11, 13, 14, 15, 16, 18, 19, 20, 21, 23, 24, 25, 27);
	/* (a0 + a3 + a5 + a6 + a7 + a10 + a12 + a13 + a14)(b0 + b3 + b5 + b6 + b7 + b10 + b12 + b13 + b14)*(X + X^2 + X^4 + X^6 + X^8 + X0 + X6 + X7 + X9 + X^21 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a3, a5, a6, a7, a10, a12, a13, a14),
				mzd_sum(t2, 9, b0, b3, b5, b6, b7, b10, b12, b13, b14), 0), X, 12,  1, 2, 4, 6, 8, 10, 16, 17, 19, 21, 23, 25);
	/* (a1 + a3 + a4 + a5 + a8 + a10 + a11 + a12)(b1 + b3 + b4 + b5 + b8 + b10 + b11 + b12)*(X^2 + X^3 + X^5 + X^7 + X^9 + X1 + X7 + X8 + X^20 + X^22 + X^24 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a3, a4, a5, a8, a10, a11, a12),
				mzd_sum(t2, 8, b1, b3, b4, b5, b8, b10, b11, b12), 0), X, 12,  2, 3, 5, 7, 9, 11, 17, 18, 20, 22, 24, 26);
	/* (a2 + a4 + a5 + a6 + a9 + a11 + a12 + a13)(b2 + b4 + b5 + b6 + b9 + b11 + b12 + b13)*(X^3 + X^4 + X^6 + X^8 + X0 + X2 + X8 + X9 + X^21 + X^23 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a4, a5, a6, a9, a11, a12, a13),
				mzd_sum(t2, 8, b2, b4, b5, b6, b9, b11, b12, b13), 0), X, 12,  3, 4, 6, 8, 10, 12, 18, 19, 21, 23, 25, 27);
	/* (a0 + a1 + a4 + a6 + a7 + a8 + a11 + a13 + a14)(b0 + b1 + b4 + b6 + b7 + b8 + b11 + b13 + b14)*(X + X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X1 + X6 + X8 + X9 + X^20 + X^21 + X^22 + X^23 + X^24 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a4, a6, a7, a8, a11, a13, a14),
				mzd_sum(t2, 9, b0, b1, b4, b6, b7, b8, b11, b13, b14), 0), X, 20,  1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26);
	/* (a0 + a2 + a3 + a4 + a7 + a9 + a10 + a11 + a14)(b0 + b2 + b3 + b4 + b7 + b9 + b10 + b11 + b14)*(X^2 + X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X1 + X2 + X7 + X9 + X^20 + X^21 + X^22 + X^23 + X^24 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a3, a4, a7, a9, a10, a11, a14),
				mzd_sum(t2, 9, b0, b2, b3, b4, b7, b9, b10, b11, b14), 0), X, 20,  2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27);
	/* (a1 + a2 + a3 + a6 + a8 + a9 + a10 + a13)(b1 + b2 + b3 + b6 + b8 + b9 + b10 + b13)*(X + X^5 + X^7 + X^9 + X1 + X2 + X6 + X^20 + X^22 + X^24 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a2, a3, a6, a8, a9, a10, a13),
				mzd_sum(t2, 8, b1, b2, b3, b6, b8, b9, b10, b13), 0), X, 12,  1, 5, 7, 9, 11, 12, 16, 20, 22, 24, 26, 27);
	/* (a0 + a3 + a4 + a5 + a7 + a10 + a11 + a12 + a14)(b0 + b3 + b4 + b5 + b7 + b10 + b11 + b12 + b14)*(X^3 + X^6 + X^9 + X0 + X1 + X2 + X8 + X^21 + X^24 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a3, a4, a5, a7, a10, a11, a12, a14),
				mzd_sum(t2, 9, b0, b3, b4, b5, b7, b10, b11, b12, b14), 0), X, 12,  3, 6, 9, 10, 11, 12, 18, 21, 24, 25, 26, 27);
	/* (a1 + a4 + a5 + a6 + a8 + a11 + a12 + a13)(b1 + b4 + b5 + b6 + b8 + b11 + b12 + b13)*(X + X^3 + X^4 + X^6 + X^7 + X^8 + X1 + X2 + X6 + X8 + X9 + X^21 + X^22 + X^23 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a4, a5, a6, a8, a11, a12, a13),
				mzd_sum(t2, 8, b1, b4, b5, b6, b8, b11, b12, b13), 0), X, 16,  1, 3, 4, 6, 7, 8, 11, 12, 16, 18, 19, 21, 22, 23, 26, 27);
	/* (a2 + a3 + a4 + a6 + a9 + a10 + a11 + a13)(b2 + b3 + b4 + b6 + b9 + b10 + b11 + b13)*(X + X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^9 + X0 + X2 + X6 + X7 + X8 + X9 + X^20 + X^21 + X^22 + X^24 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a3, a4, a6, a9, a10, a11, a13),
				mzd_sum(t2, 8, b2, b3, b4, b6, b9, b10, b11, b13), 0), X, 20,  1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 16, 17, 18, 19, 20, 21, 22, 24, 25, 27);
	/* (a0 + a1 + a3 + a6 + a7 + a8 + a10 + a13 + a14)(b0 + b1 + b3 + b6 + b7 + b8 + b10 + b13 + b14)*(X^2 + X^3 + X^5 + X^6 + X^8 + X2 + X7 + X8 + X^20 + X^21 + X^23 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a3, a6, a7, a8, a10, a13, a14),
				mzd_sum(t2, 9, b0, b1, b3, b6, b7, b8, b10, b13, b14), 0), X, 12,  2, 3, 5, 6, 8, 12, 17, 18, 20, 21, 23, 27);
	/* (a0 + a2 + a5 + a6 + a7 + a9 + a12 + a13 + a14)(b0 + b2 + b5 + b6 + b7 + b9 + b12 + b13 + b14)*(X + X^4 + X^7 + X^8 + X^9 + X0 + X6 + X9 + X^22 + X^23 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a5, a6, a7, a9, a12, a13, a14),
				mzd_sum(t2, 9, b0, b2, b5, b6, b7, b9, b12, b13, b14), 0), X, 12,  1, 4, 7, 8, 9, 10, 16, 19, 22, 23, 24, 25);
	/* (a1 + a2 + a3 + a5 + a8 + a9 + a10 + a12)(b1 + b2 + b3 + b5 + b8 + b9 + b10 + b12)*(X^2 + X^5 + X^8 + X^9 + X0 + X1 + X7 + X^20 + X^23 + X^24 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a2, a3, a5, a8, a9, a10, a12),
				mzd_sum(t2, 8, b1, b2, b3, b5, b8, b9, b10, b12), 0), X, 12,  2, 5, 8, 9, 10, 11, 17, 20, 23, 24, 25, 26);
	/* (a0 + a4 + a7 + a8 + a10 + a12 + a13 + a14)(b0 + b4 + b7 + b8 + b10 + b12 + b13 + b14)*(X + X^2 + X^4 + X^6 + X^7 + X^8 + X2 + X3 + X4 + X5 + X7 + X8 + X9 + X^20 + X^22 + X^23 + X^24 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a4, a7, a8, a10, a12, a13, a14),
				mzd_sum(t2, 8, b0, b4, b7, b8, b10, b12, b13, b14), 0), X, 18,  1, 2, 4, 6, 7, 8, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24, 27);
	/* (a1 + a4 + a5 + a7 + a9 + a10 + a11 + a12)(b1 + b4 + b5 + b7 + b9 + b10 + b11 + b12)*(X + X^2 + X^5 + X^6 + X^7 + X^9 + X0 + X4 + X5 + X9 + X^20 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a4, a5, a7, a9, a10, a11, a12),
				mzd_sum(t2, 8, b1, b4, b5, b7, b9, b10, b11, b12), 0), X, 12,  1, 2, 5, 6, 7, 9, 10, 14, 15, 19, 20, 24);
	/* (a2 + a5 + a6 + a8 + a10 + a11 + a12 + a13)(b2 + b5 + b6 + b8 + b10 + b11 + b12 + b13)*(X^2 + X^3 + X^6 + X^7 + X^8 + X0 + X1 + X5 + X6 + X^20 + X^21 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a5, a6, a8, a10, a11, a12, a13),
				mzd_sum(t2, 8, b2, b5, b6, b8, b10, b11, b12, b13), 0), X, 12,  2, 3, 6, 7, 8, 10, 11, 15, 16, 20, 21, 25);
	/* (a3 + a6 + a7 + a9 + a11 + a12 + a13 + a14)(b3 + b6 + b7 + b9 + b11 + b12 + b13 + b14)*(X^3 + X^4 + X^7 + X^8 + X^9 + X1 + X2 + X6 + X7 + X^21 + X^22 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a3, a6, a7, a9, a11, a12, a13, a14),
				mzd_sum(t2, 8, b3, b6, b7, b9, b11, b12, b13, b14), 0), X, 12,  3, 4, 7, 8, 9, 11, 12, 16, 17, 21, 22, 26);
	/* (a0 + a1 + a2 + a3 + a7 + a10 + a11 + a13)(b0 + b1 + b2 + b3 + b7 + b10 + b11 + b13)*(X^2 + X^4 + X^6 + X^9 + X0 + X2 + X5 + X7 + X^20 + X^22 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a2, a3, a7, a10, a11, a13),
				mzd_sum(t2, 8, b0, b1, b2, b3, b7, b10, b11, b13), 0), X, 12,  2, 4, 6, 9, 10, 12, 15, 17, 20, 22, 25, 26);
	/* (a0 + a1 + a5 + a8 + a9 + a11 + a13 + a14)(b0 + b1 + b5 + b8 + b9 + b11 + b13 + b14)*(X^2 + X^5 + X^6 + X^8 + X3 + X5 + X8 + X^20 + X^23 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a5, a8, a9, a11, a13, a14),
				mzd_sum(t2, 8, b0, b1, b5, b8, b9, b11, b13, b14), 0), X, 12,  2, 5, 6, 8, 13, 15, 18, 20, 23, 25, 26, 27);
	/* (a0 + a2 + a4 + a5 + a6 + a7 + a11 + a14)(b0 + b2 + b4 + b5 + b6 + b7 + b11 + b14)*(X + X^2 + X^3 + X^4 + X^5 + X^6 + X^8 + X0 + X1 + X2 + X4 + X5 + X6 + X7 + X9 + X^20 + X^21 + X^22 + X^24 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a4, a5, a6, a7, a11, a14),
				mzd_sum(t2, 8, b0, b2, b4, b5, b6, b7, b11, b14), 0), X, 20,  1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 22, 24, 26);
	/* (a1 + a3 + a4 + a5 + a6 + a10 + a13 + a14)(b1 + b3 + b4 + b5 + b6 + b10 + b13 + b14)*(X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^9 + X1 + X2 + X3 + X5 + X6 + X7 + X8 + X^20 + X^21 + X^22 + X^23 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a3, a4, a5, a6, a10, a13, a14),
				mzd_sum(t2, 8, b1, b3, b4, b5, b6, b10, b13, b14), 0), X, 20,  2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 15, 16, 17, 18, 20, 21, 22, 23, 25, 27);
	/* (a2 + a3 + a5 + a7 + a8 + a9 + a10 + a14)(b2 + b3 + b5 + b7 + b8 + b9 + b10 + b14)*(X + X^2 + X^3 + X^6 + X^9 + X1 + X3 + X4 + X5 + X6 + X8 + X9 + X^20 + X^21 + X^23 + X^24 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a3, a5, a7, a8, a9, a10, a14),
				mzd_sum(t2, 8, b2, b3, b5, b7, b8, b9, b10, b14), 0), X, 18,  1, 2, 3, 6, 9, 11, 13, 14, 15, 16, 18, 19, 20, 21, 23, 24, 26, 27);
	/* (a0 + a4 + a5 + a6 + a7 + a9 + a11 + a12)(b0 + b4 + b5 + b6 + b7 + b9 + b11 + b12)*(X^3 + X^5 + X^6 + X^7 + X0 + X1 + X2 + X3 + X5 + X7 + X9 + X^21 + X^22 + X^24 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a4, a5, a6, a7, a9, a11, a12),
				mzd_sum(t2, 8, b0, b4, b5, b6, b7, b9, b11, b12), 0), X, 16,  3, 5, 6, 7, 10, 11, 12, 13, 15, 17, 19, 21, 22, 24, 25, 26);
	/* (a1 + a5 + a6 + a7 + a8 + a10 + a12 + a13)(b1 + b5 + b6 + b7 + b8 + b10 + b12 + b13)*(X^4 + X^6 + X^7 + X^8 + X1 + X2 + X3 + X4 + X6 + X8 + X^20 + X^22 + X^23 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a5, a6, a7, a8, a10, a12, a13),
				mzd_sum(t2, 8, b1, b5, b6, b7, b8, b10, b12, b13), 0), X, 16,  4, 6, 7, 8, 11, 12, 13, 14, 16, 18, 20, 22, 23, 25, 26, 27);
	/* (a2 + a6 + a7 + a8 + a9 + a11 + a13 + a14)(b2 + b6 + b7 + b8 + b9 + b11 + b13 + b14)*(X + X^3 + X^5 + X^6 + X^7 + X^9 + X0 + X2 + X4 + X5 + X6 + X7 + X8 + X9 + X^24 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a6, a7, a8, a9, a11, a13, a14),
				mzd_sum(t2, 8, b2, b6, b7, b8, b9, b11, b13, b14), 0), X, 18,  1, 3, 5, 6, 7, 9, 10, 12, 14, 15, 16, 17, 18, 19, 24, 25, 26, 27);
	/* (a3 + a4 + a5 + a6 + a8 + a10 + a11 + a14)(b3 + b4 + b5 + b6 + b8 + b10 + b11 + b14)*(X + X^2 + X^3 + X^4 + X^7 + X1 + X5 + X7 + X9 + X^20 + X^21 + X^23 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a3, a4, a5, a6, a8, a10, a11, a14),
				mzd_sum(t2, 8, b3, b4, b5, b6, b8, b10, b11, b14), 0), X, 14,  1, 2, 3, 4, 7, 11, 15, 17, 19, 20, 21, 23, 26, 27);
	/* (a0 + a1 + a2 + a3 + a5 + a7 + a8 + a11)(b0 + b1 + b2 + b3 + b5 + b7 + b8 + b11)*(X^3 + X^4 + X^5 + X^8 + X0 + X4 + X5 + X6 + X7 + X8 + X9 + X^20 + X^21 + X^23 + X^24 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a2, a3, a5, a7, a8, a11),
				mzd_sum(t2, 8, b0, b1, b2, b3, b5, b7, b8, b11), 0), X, 16,  3, 4, 5, 8, 10, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 27);
	/* (a0 + a1 + a4 + a8 + a9 + a10 + a11 + a13)(b0 + b1 + b4 + b8 + b9 + b10 + b11 + b13)*(X + X^2 + X^4 + X^5 + X^6 + X0 + X2 + X3 + X^20 + X^22 + X^23 + X^24 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a4, a8, a9, a10, a11, a13),
				mzd_sum(t2, 8, b0, b1, b4, b8, b9, b10, b11, b13), 0), X, 14,  1, 2, 4, 5, 6, 10, 12, 13, 20, 22, 23, 24, 25, 27);
	/* (a0 + a2 + a4 + a5 + a8 + a12 + a13 + a14)(b0 + b2 + b4 + b5 + b8 + b12 + b13 + b14)*(X^2 + X^5 + X^7 + X^8 + X^9 + X0 + X3 + X^21 + X^22 + X^24 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a4, a5, a8, a12, a13, a14),
				mzd_sum(t2, 8, b0, b2, b4, b5, b8, b12, b13, b14), 0), X, 12,  2, 5, 7, 8, 9, 10, 13, 21, 22, 24, 26, 27);
	/* (a1 + a3 + a4 + a7 + a11 + a12 + a13 + a14)(b1 + b3 + b4 + b7 + b11 + b12 + b13 + b14)*(X + X^9 + X1 + X3 + X4 + X6 + X8 + X^21 + X^22 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a3, a4, a7, a11, a12, a13, a14),
				mzd_sum(t2, 8, b1, b3, b4, b7, b11, b12, b13, b14), 0), X, 10,  1, 9, 11, 13, 14, 16, 18, 21, 22, 27);
	/* (a2 + a3 + a4 + a5 + a7 + a9 + a10 + a13)(b2 + b3 + b4 + b5 + b7 + b9 + b10 + b13)*(X^2 + X^3 + X^6 + X^8 + X^9 + X1 + X2 + X5 + X7 + X9 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a3, a4, a5, a7, a9, a10, a13),
				mzd_sum(t2, 8, b2, b3, b4, b5, b7, b9, b10, b13), 0), X, 12,  2, 3, 6, 8, 9, 11, 12, 15, 17, 19, 25, 27);
	/* (a0 + a4 + a5 + a9 + a10 + a14)(b0 + b4 + b5 + b9 + b10 + b14)*(X^3 + X^4 + X^5 + X^6 + X2 + X3 + X5 + X6 + X7 + X9 + X^22 + X^23 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a4, a5, a9, a10, a14),
				mzd_sum(t2, 6, b0, b4, b5, b9, b10, b14), 0), X, 14,  3, 4, 5, 6, 12, 13, 15, 16, 17, 19, 22, 23, 25, 26);
	/* (a1 + a4 + a6 + a9 + a11 + a14)(b1 + b4 + b6 + b9 + b11 + b14)*(X^4 + X^5 + X^6 + X^7 + X3 + X4 + X6 + X7 + X8 + X^20 + X^23 + X^24 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a4, a6, a9, a11, a14),
				mzd_sum(t2, 6, b1, b4, b6, b9, b11, b14), 0), X, 14,  4, 5, 6, 7, 13, 14, 16, 17, 18, 20, 23, 24, 26, 27);
	/* (a2 + a4 + a7 + a9 + a12 + a14)(b2 + b4 + b7 + b9 + b12 + b14)*(X + X^3 + X^5 + X^7 + X0 + X3 + X4 + X5 + X6 + X7 + X9 + X^23 + X^24 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a4, a7, a9, a12, a14),
				mzd_sum(t2, 6, b2, b4, b7, b9, b12, b14), 0), X, 14,  1, 3, 5, 7, 10, 13, 14, 15, 16, 17, 19, 23, 24, 27);
	/* (a3 + a4 + a8 + a9 + a13 + a14)(b3 + b4 + b8 + b9 + b13 + b14)*(X + X^2 + X^3 + X^4 + X0 + X1 + X3 + X4 + X5 + X7 + X^20 + X^21 + X^23 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a3, a4, a8, a9, a13, a14),
				mzd_sum(t2, 6, b3, b4, b8, b9, b13, b14), 0), X, 14,  1, 2, 3, 4, 10, 11, 13, 14, 15, 17, 20, 21, 23, 24);
	/* (a0 + a1 + a2 + a3 + a5 + a6 + a7 + a8 + a10 + a11 + a12 + a13)(b0 + b1 + b2 + b3 + b5 + b6 + b7 + b8 + b10 + b11 + b12 + b13)*(X^2 + X^3 + X^4 + X^5 + X1 + X2 + X4 + X5 + X6 + X8 + X^21 + X^22 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 12, a0, a1, a2, a3, a5, a6, a7, a8, a10, a11, a12, a13),
				mzd_sum(t2, 12, b0, b1, b2, b3, b5, b6, b7, b8, b10, b11, b12, b13), 0), X, 14,  2, 3, 4, 5, 11, 12, 14, 15, 16, 18, 21, 22, 24, 25);
	/* (a0 + a1 + a5 + a6 + a10 + a11)(b0 + b1 + b5 + b6 + b10 + b11)*(X + X^2 + X^4 + X^7 + X0 + X1 + X2 + X3 + X7 + X8 + X9 + X^21 + X^22 + X^23 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a5, a6, a10, a11),
				mzd_sum(t2, 6, b0, b1, b5, b6, b10, b11), 0), X, 16,  1, 2, 4, 7, 10, 11, 12, 13, 17, 18, 19, 21, 22, 23, 25, 27);
	/* (a0 + a2 + a5 + a7 + a10 + a12)(b0 + b2 + b5 + b7 + b10 + b12)*(X + X^5 + X0 + X2 + X3 + X6 + X7 + X8 + X^20 + X^22 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a2, a5, a7, a10, a12),
				mzd_sum(t2, 6, b0, b2, b5, b7, b10, b12), 0), X, 12,  1, 5, 10, 12, 13, 16, 17, 18, 20, 22, 23, 25);
	/* (a1 + a3 + a6 + a8 + a11 + a13)(b1 + b3 + b6 + b8 + b11 + b13)*(X^2 + X^6 + X1 + X3 + X4 + X7 + X8 + X9 + X^21 + X^23 + X^24 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a6, a8, a11, a13),
				mzd_sum(t2, 6, b1, b3, b6, b8, b11, b13), 0), X, 12,  2, 6, 11, 13, 14, 17, 18, 19, 21, 23, 24, 26);
	/* (a2 + a3 + a7 + a8 + a12 + a13)(b2 + b3 + b7 + b8 + b12 + b13)*(X^2 + X^3 + X^6 + X^7 + X1 + X2 + X3 + X5 + X7 + X^20 + X^21 + X^22 + X^23 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a3, a7, a8, a12, a13),
				mzd_sum(t2, 6, b2, b3, b7, b8, b12, b13), 0), X, 16,  2, 3, 6, 7, 11, 12, 13, 15, 17, 20, 21, 22, 23, 25, 26, 27);
	/* (a0 + a5 + a8 + a10 + a11 + a14)(b0 + b5 + b8 + b10 + b11 + b14)*(X + X^3 + X^8 + X0 + X6 + X8 + X^23 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a0, a5, a8, a10, a11, a14),
				mzd_sum(t2, 6, b0, b5, b8, b10, b11, b14), 0), X, 8,  1, 3, 8, 10, 16, 18, 23, 25);
	/* (a1 + a6 + a9 + a11 + a12)(b1 + b6 + b9 + b11 + b12)*(X + X^3 + X^4 + X^5 + X^8 + X0 + X1 + X2 + X6 + X8 + X9 + X^20 + X^23 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 5, a1, a6, a9, a11, a12),
				mzd_sum(t2, 5, b1, b6, b9, b11, b12), 0), X, 16,  1, 3, 4, 5, 8, 10, 11, 12, 16, 18, 19, 20, 23, 25, 26, 27);
	/* (a3 + a6 + a8 + a9 + a12 + a13 + a14)(b3 + b6 + b8 + b9 + b12 + b13 + b14)*(X + X^2 + X^5 + X^8 + X^9 + X2 + X6 + X7 + X^20 + X^23 + X^24 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a3, a6, a8, a9, a12, a13, a14),
				mzd_sum(t2, 7, b3, b6, b8, b9, b12, b13, b14), 0), X, 12,  1, 2, 5, 8, 9, 12, 16, 17, 20, 23, 24, 27);
	/* (a4 + a7 + a9 + a10 + a13 + a14)(b4 + b7 + b9 + b10 + b13 + b14)*(X^3 + X^5 + X0 + X2 + X8 + X^20 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a4, a7, a9, a10, a13, a14),
				mzd_sum(t2, 6, b4, b7, b9, b10, b13, b14), 0), X, 8,  3, 5, 10, 12, 18, 20, 25, 27);
	/* (a0 + a4 + a5 + a7 + a8 + a9 + a11 + a13)(b0 + b4 + b5 + b7 + b8 + b9 + b11 + b13)*(X^3 + X0 + X8 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a4, a5, a7, a8, a9, a11, a13),
				mzd_sum(t2, 8, b0, b4, b5, b7, b8, b9, b11, b13), 0), X, 4,  3, 10, 18, 25);
	/* (a0 + a1 + a5 + a6 + a8 + a9 + a10 + a12 + a14)(b0 + b1 + b5 + b6 + b8 + b9 + b10 + b12 + b14)*(X^2 + X^9 + X7 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a5, a6, a8, a9, a10, a12, a14),
				mzd_sum(t2, 9, b0, b1, b5, b6, b8, b9, b10, b12, b14), 0), X, 4,  2, 9, 17, 24);
	/* (a3 + a4 + a6 + a7 + a8 + a10 + a12)(b3 + b4 + b6 + b7 + b8 + b10 + b12)*(X^4 + X1 + X9 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a3, a4, a6, a7, a8, a10, a12),
				mzd_sum(t2, 7, b3, b4, b6, b7, b8, b10, b12), 0), X, 4,  4, 11, 19, 26);
	/* (a1 + a2 + a4 + a5 + a6 + a8 + a10)(b1 + b2 + b4 + b5 + b6 + b8 + b10)*(X + X^3 + X^5 + X^8 + X0 + X2 + X6 + X8 + X^20 + X^23 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a1, a2, a4, a5, a6, a8, a10),
				mzd_sum(t2, 7, b1, b2, b4, b5, b6, b8, b10), 0), X, 12,  1, 3, 5, 8, 10, 12, 16, 18, 20, 23, 25, 27);
	/* (a0 + a2 + a3 + a6 + a7 + a8 + a9 + a10 + a14)(b0 + b2 + b3 + b6 + b7 + b8 + b9 + b10 + b14)*(X^2 + X^3 + X^4 + X^5 + X^9 + X0 + X1 + X2 + X7 + X8 + X9 + X^20 + X^24 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a3, a6, a7, a8, a9, a10, a14),
				mzd_sum(t2, 9, b0, b2, b3, b6, b7, b8, b9, b10, b14), 0), X, 16,  2, 3, 4, 5, 9, 10, 11, 12, 17, 18, 19, 20, 24, 25, 26, 27);
	/* (a0 + a1 + a3 + a4 + a5 + a7 + a9 + a14)(b0 + b1 + b3 + b4 + b5 + b7 + b9 + b14)*(X + X^2 + X^3 + X^4 + X^8 + X^9 + X0 + X1 + X6 + X7 + X8 + X9 + X^23 + X^24 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a3, a4, a5, a7, a9, a14),
				mzd_sum(t2, 8, b0, b1, b3, b4, b5, b7, b9, b14), 0), X, 16,  1, 2, 3, 4, 8, 9, 10, 11, 16, 17, 18, 19, 23, 24, 25, 26);
	/* (a0 + a1 + a2 + a4 + a6 + a11 + a14)(b0 + b1 + b2 + b4 + b6 + b11 + b14)*(X + X^2 + X^3 + X^4 + X^5 + X^8 + X^9 + X0 + X1 + X2 + X6 + X7 + X8 + X9 + X^20 + X^23 + X^24 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a2, a4, a6, a11, a14),
				mzd_sum(t2, 7, b0, b1, b2, b4, b6, b11, b14), 0), X, 20,  1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27);
	/* (a0 + a2 + a3 + a4 + a6 + a8 + a13)(b0 + b2 + b3 + b4 + b6 + b8 + b13)*(X + X^2 + X^4 + X^5 + X^8 + X^9 + X1 + X2 + X6 + X7 + X9 + X^20 + X^23 + X^24 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a3, a4, a6, a8, a13),
				mzd_sum(t2, 7, b0, b2, b3, b4, b6, b8, b13), 0), X, 16,  1, 2, 4, 5, 8, 9, 11, 12, 16, 17, 19, 20, 23, 24, 26, 27);
	/* (a0 + a1 + a2 + a3 + a4 + a8 + a9 + a11 + a12 + a13)(b0 + b1 + b2 + b3 + b4 + b8 + b9 + b11 + b12 + b13)*(X + X^4 + X^8 + X1 + X6 + X9 + X^23 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a1, a2, a3, a4, a8, a9, a11, a12, a13),
				mzd_sum(t2, 10, b0, b1, b2, b3, b4, b8, b9, b11, b12, b13), 0), X, 8,  1, 4, 8, 11, 16, 19, 23, 26);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba16(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
	/* using 64 multiplications and three temporaries (not optimal!) */
	const mzd_t *a0  = A[ 0];
	const mzd_t *a1  = A[ 1];
	const mzd_t *a2  = A[ 2];
	const mzd_t *a3  = A[ 3];
	const mzd_t *a4  = A[ 4];
	const mzd_t *a5  = A[ 5];
	const mzd_t *a6  = A[ 6];
	const mzd_t *a7  = A[ 7];
	const mzd_t *a8  = A[ 8];
	const mzd_t *a9  = A[ 9];
	const mzd_t *a10 = A[10];
	const mzd_t *a11 = A[11];
	const mzd_t *a12 = A[12];
	const mzd_t *a13 = A[13];
	const mzd_t *a14 = A[14];
	const mzd_t *a15 = A[15];

	const mzd_t *b0  = B[ 0];
	const mzd_t *b1  = B[ 1];
	const mzd_t *b2  = B[ 2];
	const mzd_t *b3  = B[ 3];
	const mzd_t *b4  = B[ 4];
	const mzd_t *b5  = B[ 5];
	const mzd_t *b6  = B[ 6];
	const mzd_t *b7  = B[ 7];
	const mzd_t *b8  = B[ 8];
	const mzd_t *b9  = B[ 9];
	const mzd_t *b10 = B[10];
	const mzd_t *b11 = B[11];
	const mzd_t *b12 = B[12];
	const mzd_t *b13 = B[13];
	const mzd_t *b14 = B[14];
	const mzd_t *b15 = B[15];

	mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
	mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
	mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

	/* (a0 + a2 + a3 + a5 + a6 + a8 + a9 + a11 + a12 + a14 + a15)(b0 + b2 + b3 + b5 + b6 + b8 + b9 + b11 + b12 + b14 + b15)*(X^2 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X1 + X2 + X5 + X6 + X7 + X8 + X9 + X^20 + X^23 + X^26 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 11, a0, a2, a3, a5, a6, a8, a9, a11, a12, a14, a15),
				mzd_sum(t2, 11, b0, b2, b3, b5, b6, b8, b9, b11, b12, b14, b15), 0), X, 18,  2, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19, 20, 23, 26, 28);
	/* (a1 + a2 + a4 + a5 + a7 + a8 + a10 + a11 + a13 + a14)(b1 + b2 + b4 + b5 + b7 + b8 + b10 + b11 + b13 + b14)*(X^2 + X^3 + X^4 + X^6 + X^8 + X0 + X2 + X3 + X4 + X6 + X8 + X^20 + X^21 + X^22 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a1, a2, a4, a5, a7, a8, a10, a11, a13, a14),
				mzd_sum(t2, 10, b1, b2, b4, b5, b7, b8, b10, b11, b13, b14), 0), X, 16,  2, 3, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 21, 22, 26, 27);
	/* (a0 + a1 + a3 + a4 + a6 + a7 + a9 + a10 + a12 + a13 + a15)(b0 + b1 + b3 + b4 + b6 + b7 + b9 + b10 + b12 + b13 + b15)*(X^3 + X^4 + X^5 + X^7 + X^9 + X1 + X3 + X4 + X5 + X7 + X9 + X^21 + X^22 + X^23 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 11, a0, a1, a3, a4, a6, a7, a9, a10, a12, a13, a15),
				mzd_sum(t2, 11, b0, b1, b3, b4, b6, b7, b9, b10, b12, b13, b15), 0), X, 16,  3, 4, 5, 7, 9, 11, 13, 14, 15, 17, 19, 21, 22, 23, 27, 28);
	/* (a0)(b0)*(1 + X + X^2 + X^3 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X2 + X3 + X5 + X6 + X7 + X8 + X^20 + X^21 + X^22 + X^23 + X^24 + X^25 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a0, b0, 0), X, 24,  0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25, 27, 28);
	/* (a1)(b1)*(X + X^3 + X^6 + X^8 + X0 + X3 + X6 + X8 + X^21 + X^23 + X^25 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a1, b1, 0), X, 12,  1, 3, 6, 8, 10, 13, 16, 18, 21, 23, 25, 28);
	/* (a0 + a1)(b0 + b1)*(X + X^3 + X^6 + X^8 + X0 + X3 + X6 + X8 + X^21 + X^23 + X^25 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
				mzd_sum(t2, 2, b0, b1), 0), X, 12,  1, 3, 6, 8, 10, 13, 16, 18, 21, 23, 25, 28);
	/* (a0 + a2 + a4 + a6 + a8 + a10 + a12 + a14)(b0 + b2 + b4 + b6 + b8 + b10 + b12 + b14)*(X^2 + X^3 + X^7 + X^8 + X1 + X2 + X3 + X7 + X8 + X^22 + X^23 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a2, a4, a6, a8, a10, a12, a14),
				mzd_sum(t2, 8, b0, b2, b4, b6, b8, b10, b12, b14), 0), X, 14,  2, 3, 7, 8, 11, 12, 13, 17, 18, 22, 23, 26, 27, 28);
	/* (a1 + a3 + a5 + a7 + a9 + a11 + a13 + a15)(b1 + b3 + b5 + b7 + b9 + b11 + b13 + b15)*(X^2 + X^3 + X^7 + X^8 + X1 + X2 + X3 + X7 + X8 + X^22 + X^23 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a3, a5, a7, a9, a11, a13, a15),
				mzd_sum(t2, 8, b1, b3, b5, b7, b9, b11, b13, b15), 0), X, 14,  2, 3, 7, 8, 11, 12, 13, 17, 18, 22, 23, 26, 27, 28);
	/* (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15)(b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10 + b11 + b12 + b13 + b14 + b15)*(X^2 + X^7 + X1 + X3 + X4 + X5 + X6 + X8 + X9 + X^20 + X^21 + X^23 + X^24 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 16, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15),
				mzd_sum(t2, 16, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15), 0), X, 15,  2, 7, 11, 13, 14, 15, 16, 18, 19, 20, 21, 23, 24, 25, 27);
	/* (a15)(b15)*(X^2 + X^3 + X^4 + X^5 + X^7 + X^8 + X^9 + X0 + X1 + X2 + X4 + X5 + X7 + X8 + X9 + X^20 + X^22 + X^23 + X^24 + X^25 + X^26 + X^27 + X^29 + X^30) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a15, b15, 0), X, 24,  2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30);
	/* (a14)(b14)*(X^2 + X^4 + X^7 + X^9 + X1 + X4 + X7 + X9 + X^22 + X^24 + X^26 + X^29) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, a14, b14, 0), X, 12,  2, 4, 7, 9, 11, 14, 17, 19, 22, 24, 26, 29);
	/* (a14 + a15)(b14 + b15)*(X^2 + X^4 + X^7 + X^9 + X1 + X4 + X7 + X9 + X^22 + X^24 + X^26 + X^29) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a14, a15),
				mzd_sum(t2, 2, b14, b15), 0), X, 12,  2, 4, 7, 9, 11, 14, 17, 19, 22, 24, 26, 29);
	/* (a0 + a3 + a5 + a6 + a7 + a10 + a12 + a13 + a14)(b0 + b3 + b5 + b6 + b7 + b10 + b12 + b13 + b14)*(X^2 + X^3 + X^4 + X3 + X7 + X8 + X9 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a3, a5, a6, a7, a10, a12, a13, a14),
				mzd_sum(t2, 9, b0, b3, b5, b6, b7, b10, b12, b13, b14), 0), X, 8,  2, 3, 4, 13, 17, 18, 19, 28);
	/* (a1 + a3 + a4 + a5 + a8 + a10 + a11 + a12 + a15)(b1 + b3 + b4 + b5 + b8 + b10 + b11 + b12 + b15)*(X^2 + X^3 + X^5 + X^7 + X^9 + X1 + X7 + X8 + X^20 + X^22 + X^24 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a1, a3, a4, a5, a8, a10, a11, a12, a15),
				mzd_sum(t2, 9, b1, b3, b4, b5, b8, b10, b11, b12, b15), 0), X, 12,  2, 3, 5, 7, 9, 11, 17, 18, 20, 22, 24, 26);
	/* (a2 + a4 + a5 + a6 + a9 + a11 + a12 + a13)(b2 + b4 + b5 + b6 + b9 + b11 + b12 + b13)*(X^3 + X^4 + X^6 + X^8 + X0 + X2 + X8 + X9 + X^21 + X^23 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a4, a5, a6, a9, a11, a12, a13),
				mzd_sum(t2, 8, b2, b4, b5, b6, b9, b11, b12, b13), 0), X, 12,  3, 4, 6, 8, 10, 12, 18, 19, 21, 23, 25, 27);
	/* (a0 + a1 + a4 + a6 + a7 + a8 + a11 + a13 + a14 + a15)(b0 + b1 + b4 + b6 + b7 + b8 + b11 + b13 + b14 + b15)*(X^4 + X^5 + X^7 + X^9 + X1 + X3 + X9 + X^20 + X^22 + X^24 + X^26 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a1, a4, a6, a7, a8, a11, a13, a14, a15),
				mzd_sum(t2, 10, b0, b1, b4, b6, b7, b8, b11, b13, b14, b15), 0), X, 12,  4, 5, 7, 9, 11, 13, 19, 20, 22, 24, 26, 28);
	/* (a0 + a2 + a3 + a4 + a7 + a9 + a10 + a11 + a14)(b0 + b2 + b3 + b4 + b7 + b9 + b10 + b11 + b14)*(X^2 + X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X1 + X2 + X7 + X9 + X^20 + X^21 + X^22 + X^23 + X^24 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a3, a4, a7, a9, a10, a11, a14),
				mzd_sum(t2, 9, b0, b2, b3, b4, b7, b9, b10, b11, b14), 0), X, 20,  2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27);
	/* (a1 + a2 + a3 + a6 + a8 + a9 + a10 + a13 + a15)(b1 + b2 + b3 + b6 + b8 + b9 + b10 + b13 + b15)*(X^3 + X^5 + X^6 + X^7 + X^8 + X^9 + X0 + X1 + X2 + X3 + X8 + X^20 + X^21 + X^22 + X^23 + X^24 + X^25 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a1, a2, a3, a6, a8, a9, a10, a13, a15),
				mzd_sum(t2, 9, b1, b2, b3, b6, b8, b9, b10, b13, b15), 0), X, 20,  3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28);
	/* (a0 + a3 + a4 + a5 + a7 + a10 + a11 + a12 + a14)(b0 + b3 + b4 + b5 + b7 + b10 + b11 + b12 + b14)*(X^3 + X^6 + X^9 + X0 + X1 + X2 + X8 + X^21 + X^24 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a3, a4, a5, a7, a10, a11, a12, a14),
				mzd_sum(t2, 9, b0, b3, b4, b5, b7, b10, b11, b12, b14), 0), X, 12,  3, 6, 9, 10, 11, 12, 18, 21, 24, 25, 26, 27);
	/* (a1 + a4 + a5 + a6 + a8 + a11 + a12 + a13 + a15)(b1 + b4 + b5 + b6 + b8 + b11 + b12 + b13 + b15)*(X^4 + X^7 + X0 + X1 + X2 + X3 + X9 + X^22 + X^25 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a1, a4, a5, a6, a8, a11, a12, a13, a15),
				mzd_sum(t2, 9, b1, b4, b5, b6, b8, b11, b12, b13, b15), 0), X, 12,  4, 7, 10, 11, 12, 13, 19, 22, 25, 26, 27, 28);
	/* (a2 + a3 + a4 + a6 + a9 + a10 + a11 + a13)(b2 + b3 + b4 + b6 + b9 + b10 + b11 + b13)*(X^2 + X^4 + X^5 + X^7 + X^8 + X^9 + X2 + X3 + X7 + X9 + X^20 + X^22 + X^23 + X^24 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a3, a4, a6, a9, a10, a11, a13),
				mzd_sum(t2, 8, b2, b3, b4, b6, b9, b10, b11, b13), 0), X, 16,  2, 4, 5, 7, 8, 9, 12, 13, 17, 19, 20, 22, 23, 24, 27, 28);
	/* (a0 + a1 + a3 + a6 + a7 + a8 + a10 + a13 + a14 + a15)(b0 + b1 + b3 + b6 + b7 + b8 + b10 + b13 + b14 + b15)*(X^2 + X^3 + X^5 + X^6 + X^8 + X2 + X7 + X8 + X^20 + X^21 + X^23 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a1, a3, a6, a7, a8, a10, a13, a14, a15),
				mzd_sum(t2, 10, b0, b1, b3, b6, b7, b8, b10, b13, b14, b15), 0), X, 12,  2, 3, 5, 6, 8, 12, 17, 18, 20, 21, 23, 27);
	/* (a0 + a2 + a5 + a6 + a7 + a9 + a12 + a13 + a14)(b0 + b2 + b5 + b6 + b7 + b9 + b12 + b13 + b14)*(X^3 + X^4 + X^6 + X^7 + X^9 + X3 + X8 + X9 + X^21 + X^22 + X^24 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a5, a6, a7, a9, a12, a13, a14),
				mzd_sum(t2, 9, b0, b2, b5, b6, b7, b9, b12, b13, b14), 0), X, 12,  3, 4, 6, 7, 9, 13, 18, 19, 21, 22, 24, 28);
	/* (a1 + a2 + a3 + a5 + a8 + a9 + a10 + a12 + a15)(b1 + b2 + b3 + b5 + b8 + b9 + b10 + b12 + b15)*(X^2 + X^5 + X^8 + X^9 + X0 + X1 + X7 + X^20 + X^23 + X^24 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a1, a2, a3, a5, a8, a9, a10, a12, a15),
				mzd_sum(t2, 9, b1, b2, b3, b5, b8, b9, b10, b12, b15), 0), X, 12,  2, 5, 8, 9, 10, 11, 17, 20, 23, 24, 25, 26);
	/* (a0 + a4 + a7 + a8 + a10 + a12 + a13 + a14 + a15)(b0 + b4 + b7 + b8 + b10 + b12 + b13 + b14 + b15)*(X^2 + X^3 + X^4 + X^7 + X0 + X2 + X4 + X5 + X6 + X7 + X9 + X^20 + X^21 + X^22 + X^24 + X^25 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a4, a7, a8, a10, a12, a13, a14, a15),
				mzd_sum(t2, 9, b0, b4, b7, b8, b10, b12, b13, b14, b15), 0), X, 18,  2, 3, 4, 7, 10, 12, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 27, 28);
	/* (a1 + a4 + a5 + a7 + a9 + a10 + a11 + a12)(b1 + b4 + b5 + b7 + b9 + b10 + b11 + b12)*(X^2 + X^3 + X^5 + X^7 + X^8 + X^9 + X3 + X4 + X5 + X6 + X8 + X9 + X^20 + X^21 + X^23 + X^24 + X^25 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a4, a5, a7, a9, a10, a11, a12),
				mzd_sum(t2, 8, b1, b4, b5, b7, b9, b10, b11, b12), 0), X, 18,  2, 3, 5, 7, 8, 9, 13, 14, 15, 16, 18, 19, 20, 21, 23, 24, 25, 28);
	/* (a2 + a5 + a6 + a8 + a10 + a11 + a12 + a13)(b2 + b5 + b6 + b8 + b10 + b11 + b12 + b13)*(X^2 + X^3 + X^6 + X^7 + X^8 + X0 + X1 + X5 + X6 + X^20 + X^21 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a5, a6, a8, a10, a11, a12, a13),
				mzd_sum(t2, 8, b2, b5, b6, b8, b10, b11, b12, b13), 0), X, 12,  2, 3, 6, 7, 8, 10, 11, 15, 16, 20, 21, 25);
	/* (a3 + a6 + a7 + a9 + a11 + a12 + a13 + a14)(b3 + b6 + b7 + b9 + b11 + b12 + b13 + b14)*(X^3 + X^4 + X^7 + X^8 + X^9 + X1 + X2 + X6 + X7 + X^21 + X^22 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a3, a6, a7, a9, a11, a12, a13, a14),
				mzd_sum(t2, 8, b3, b6, b7, b9, b11, b12, b13, b14), 0), X, 12,  3, 4, 7, 8, 9, 11, 12, 16, 17, 21, 22, 26);
	/* (a0 + a1 + a2 + a3 + a7 + a10 + a11 + a13 + a15)(b0 + b1 + b2 + b3 + b7 + b10 + b11 + b13 + b15)*(X^2 + X^4 + X^6 + X^9 + X0 + X2 + X5 + X7 + X^20 + X^22 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a2, a3, a7, a10, a11, a13, a15),
				mzd_sum(t2, 9, b0, b1, b2, b3, b7, b10, b11, b13, b15), 0), X, 12,  2, 4, 6, 9, 10, 12, 15, 17, 20, 22, 25, 26);
	/* (a0 + a1 + a5 + a8 + a9 + a11 + a13 + a14 + a15)(b0 + b1 + b5 + b8 + b9 + b11 + b13 + b14 + b15)*(X^2 + X^5 + X^6 + X^8 + X3 + X5 + X8 + X^20 + X^23 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a5, a8, a9, a11, a13, a14, a15),
				mzd_sum(t2, 9, b0, b1, b5, b8, b9, b11, b13, b14, b15), 0), X, 12,  2, 5, 6, 8, 13, 15, 18, 20, 23, 25, 26, 27);
	/* (a0 + a2 + a4 + a5 + a6 + a7 + a11 + a14 + a15)(b0 + b2 + b4 + b5 + b6 + b7 + b11 + b14 + b15)*(X^2 + X^4 + X^5 + X1 + X2 + X3 + X4 + X5 + X7 + X8 + X9 + X^20 + X^22 + X^23 + X^24 + X^25 + X^26 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a4, a5, a6, a7, a11, a14, a15),
				mzd_sum(t2, 9, b0, b2, b4, b5, b6, b7, b11, b14, b15), 0), X, 18,  2, 4, 5, 11, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 28);
	/* (a1 + a3 + a4 + a5 + a6 + a10 + a13 + a14)(b1 + b3 + b4 + b5 + b6 + b10 + b13 + b14)*(X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^9 + X1 + X2 + X3 + X5 + X6 + X7 + X8 + X^20 + X^21 + X^22 + X^23 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a3, a4, a5, a6, a10, a13, a14),
				mzd_sum(t2, 8, b1, b3, b4, b5, b6, b10, b13, b14), 0), X, 20,  2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 15, 16, 17, 18, 20, 21, 22, 23, 25, 27);
	/* (a2 + a3 + a5 + a7 + a8 + a9 + a10 + a14)(b2 + b3 + b5 + b7 + b8 + b9 + b10 + b14)*(X^2 + X^8 + X^9 + X0 + X1 + X4 + X5 + X9 + X^20 + X^24 + X^25 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a3, a5, a7, a8, a9, a10, a14),
				mzd_sum(t2, 8, b2, b3, b5, b7, b8, b9, b10, b14), 0), X, 14,  2, 8, 9, 10, 11, 14, 15, 19, 20, 24, 25, 26, 27, 28);
	/* (a0 + a4 + a5 + a6 + a7 + a9 + a11 + a12 + a15)(b0 + b4 + b5 + b6 + b7 + b9 + b11 + b12 + b15)*(X^3 + X^5 + X^6 + X^7 + X0 + X1 + X2 + X3 + X5 + X7 + X9 + X^21 + X^22 + X^24 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a4, a5, a6, a7, a9, a11, a12, a15),
				mzd_sum(t2, 9, b0, b4, b5, b6, b7, b9, b11, b12, b15), 0), X, 16,  3, 5, 6, 7, 10, 11, 12, 13, 15, 17, 19, 21, 22, 24, 25, 26);
	/* (a1 + a5 + a6 + a7 + a8 + a10 + a12 + a13)(b1 + b5 + b6 + b7 + b8 + b10 + b12 + b13)*(X^4 + X^6 + X^7 + X^8 + X1 + X2 + X3 + X4 + X6 + X8 + X^20 + X^22 + X^23 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a5, a6, a7, a8, a10, a12, a13),
				mzd_sum(t2, 8, b1, b5, b6, b7, b8, b10, b12, b13), 0), X, 16,  4, 6, 7, 8, 11, 12, 13, 14, 16, 18, 20, 22, 23, 25, 26, 27);
	/* (a2 + a6 + a7 + a8 + a9 + a11 + a13 + a14)(b2 + b6 + b7 + b8 + b9 + b11 + b13 + b14)*(X^5 + X^7 + X^8 + X^9 + X2 + X3 + X4 + X5 + X7 + X9 + X^21 + X^23 + X^24 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a6, a7, a8, a9, a11, a13, a14),
				mzd_sum(t2, 8, b2, b6, b7, b8, b9, b11, b13, b14), 0), X, 16,  5, 7, 8, 9, 12, 13, 14, 15, 17, 19, 21, 23, 24, 26, 27, 28);
	/* (a3 + a4 + a5 + a6 + a8 + a10 + a11 + a14)(b3 + b4 + b5 + b6 + b8 + b10 + b11 + b14)*(X^2 + X^4 + X^6 + X^7 + X^8 + X0 + X1 + X3 + X5 + X6 + X7 + X8 + X9 + X^20 + X^25 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a3, a4, a5, a6, a8, a10, a11, a14),
				mzd_sum(t2, 8, b3, b4, b5, b6, b8, b10, b11, b14), 0), X, 18,  2, 4, 6, 7, 8, 10, 11, 13, 15, 16, 17, 18, 19, 20, 25, 26, 27, 28);
	/* (a0 + a1 + a2 + a3 + a5 + a7 + a8 + a11 + a15)(b0 + b1 + b2 + b3 + b5 + b7 + b8 + b11 + b15)*(X^3 + X^4 + X^5 + X^8 + X0 + X4 + X5 + X6 + X7 + X8 + X9 + X^20 + X^21 + X^23 + X^24 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a2, a3, a5, a7, a8, a11, a15),
				mzd_sum(t2, 9, b0, b1, b2, b3, b5, b7, b8, b11, b15), 0), X, 16,  3, 4, 5, 8, 10, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 27);
	/* (a0 + a1 + a4 + a8 + a9 + a10 + a11 + a13 + a15)(b0 + b1 + b4 + b8 + b9 + b10 + b11 + b13 + b15)*(X^2 + X^3 + X^4 + X^5 + X^8 + X2 + X6 + X8 + X^20 + X^21 + X^22 + X^24 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a4, a8, a9, a10, a11, a13, a15),
				mzd_sum(t2, 9, b0, b1, b4, b8, b9, b10, b11, b13, b15), 0), X, 14,  2, 3, 4, 5, 8, 12, 16, 18, 20, 21, 22, 24, 27, 28);
	/* (a0 + a2 + a4 + a5 + a8 + a12 + a13 + a14 + a15)(b0 + b2 + b4 + b5 + b8 + b12 + b13 + b14 + b15)*(X^2 + X^5 + X^7 + X^8 + X^9 + X0 + X3 + X^21 + X^22 + X^24 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a2, a4, a5, a8, a12, a13, a14, a15),
				mzd_sum(t2, 9, b0, b2, b4, b5, b8, b12, b13, b14, b15), 0), X, 12,  2, 5, 7, 8, 9, 10, 13, 21, 22, 24, 26, 27);
	/* (a1 + a3 + a4 + a7 + a11 + a12 + a13 + a14)(b1 + b3 + b4 + b7 + b11 + b12 + b13 + b14)*(X^3 + X^6 + X^8 + X^9 + X0 + X1 + X4 + X^22 + X^23 + X^25 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a3, a4, a7, a11, a12, a13, a14),
				mzd_sum(t2, 8, b1, b3, b4, b7, b11, b12, b13, b14), 0), X, 12,  3, 6, 8, 9, 10, 11, 14, 22, 23, 25, 27, 28);
	/* (a2 + a3 + a4 + a5 + a7 + a9 + a10 + a13)(b2 + b3 + b4 + b5 + b7 + b9 + b10 + b13)*(X^2 + X^3 + X^6 + X^8 + X^9 + X1 + X2 + X5 + X7 + X9 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a2, a3, a4, a5, a7, a9, a10, a13),
				mzd_sum(t2, 8, b2, b3, b4, b5, b7, b9, b10, b13), 0), X, 12,  2, 3, 6, 8, 9, 11, 12, 15, 17, 19, 25, 27);
	/* (a0 + a4 + a5 + a9 + a10 + a14 + a15)(b0 + b4 + b5 + b9 + b10 + b14 + b15)*(X^3 + X^4 + X^5 + X^6 + X2 + X3 + X5 + X6 + X7 + X9 + X^22 + X^23 + X^25 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a4, a5, a9, a10, a14, a15),
				mzd_sum(t2, 7, b0, b4, b5, b9, b10, b14, b15), 0), X, 14,  3, 4, 5, 6, 12, 13, 15, 16, 17, 19, 22, 23, 25, 26);
	/* (a1 + a4 + a6 + a9 + a11 + a14)(b1 + b4 + b6 + b9 + b11 + b14)*(X^4 + X^5 + X^6 + X^7 + X3 + X4 + X6 + X7 + X8 + X^20 + X^23 + X^24 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a4, a6, a9, a11, a14),
				mzd_sum(t2, 6, b1, b4, b6, b9, b11, b14), 0), X, 14,  4, 5, 6, 7, 13, 14, 16, 17, 18, 20, 23, 24, 26, 27);
	/* (a2 + a4 + a7 + a9 + a12 + a14)(b2 + b4 + b7 + b9 + b12 + b14)*(X^5 + X^6 + X^7 + X^8 + X4 + X5 + X7 + X8 + X9 + X^21 + X^24 + X^25 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a4, a7, a9, a12, a14),
				mzd_sum(t2, 6, b2, b4, b7, b9, b12, b14), 0), X, 14,  5, 6, 7, 8, 14, 15, 17, 18, 19, 21, 24, 25, 27, 28);
	/* (a3 + a4 + a8 + a9 + a13 + a14)(b3 + b4 + b8 + b9 + b13 + b14)*(X^2 + X^4 + X^6 + X^8 + X1 + X4 + X5 + X6 + X7 + X8 + X^20 + X^24 + X^25 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a3, a4, a8, a9, a13, a14),
				mzd_sum(t2, 6, b3, b4, b8, b9, b13, b14), 0), X, 14,  2, 4, 6, 8, 11, 14, 15, 16, 17, 18, 20, 24, 25, 28);
	/* (a0 + a1 + a2 + a3 + a5 + a6 + a7 + a8 + a10 + a11 + a12 + a13 + a15)(b0 + b1 + b2 + b3 + b5 + b6 + b7 + b8 + b10 + b11 + b12 + b13 + b15)*(X^2 + X^3 + X^4 + X^5 + X1 + X2 + X4 + X5 + X6 + X8 + X^21 + X^22 + X^24 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 13, a0, a1, a2, a3, a5, a6, a7, a8, a10, a11, a12, a13, a15),
				mzd_sum(t2, 13, b0, b1, b2, b3, b5, b6, b7, b8, b10, b11, b12, b13, b15), 0), X, 14,  2, 3, 4, 5, 11, 12, 14, 15, 16, 18, 21, 22, 24, 25);
	/* (a0 + a1 + a5 + a6 + a10 + a11 + a15)(b0 + b1 + b5 + b6 + b10 + b11 + b15)*(X^2 + X^3 + X^4 + X^6 + X^7 + X^8 + X1 + X2 + X6 + X7 + X9 + X^22 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a5, a6, a10, a11, a15),
				mzd_sum(t2, 7, b0, b1, b5, b6, b10, b11, b15), 0), X, 14,  2, 3, 4, 6, 7, 8, 11, 12, 16, 17, 19, 22, 27, 28);
	/* (a0 + a2 + a5 + a7 + a10 + a12 + a15)(b0 + b2 + b5 + b7 + b10 + b12 + b15)*(X^3 + X^5 + X^6 + X^8 + X2 + X7 + X^20 + X^21 + X^22 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a5, a7, a10, a12, a15),
				mzd_sum(t2, 7, b0, b2, b5, b7, b10, b12, b15), 0), X, 10,  3, 5, 6, 8, 12, 17, 20, 21, 22, 28);
	/* (a1 + a3 + a6 + a8 + a11 + a13)(b1 + b3 + b6 + b8 + b11 + b13)*(X^2 + X^6 + X1 + X3 + X4 + X7 + X8 + X9 + X^21 + X^23 + X^24 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a3, a6, a8, a11, a13),
				mzd_sum(t2, 6, b1, b3, b6, b8, b11, b13), 0), X, 12,  2, 6, 11, 13, 14, 17, 18, 19, 21, 23, 24, 26);
	/* (a2 + a3 + a7 + a8 + a12 + a13)(b2 + b3 + b7 + b8 + b12 + b13)*(X^2 + X^3 + X^6 + X^7 + X1 + X2 + X3 + X5 + X7 + X^20 + X^21 + X^22 + X^23 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a2, a3, a7, a8, a12, a13),
				mzd_sum(t2, 6, b2, b3, b7, b8, b12, b13), 0), X, 16,  2, 3, 6, 7, 11, 12, 13, 15, 17, 20, 21, 22, 23, 25, 26, 27);
	/* (a0 + a5 + a8 + a10 + a11 + a14 + a15)(b0 + b5 + b8 + b10 + b11 + b14 + b15)*(X^6 + X3 + X^21 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a5, a8, a10, a11, a14, a15),
				mzd_sum(t2, 7, b0, b5, b8, b10, b11, b14, b15), 0), X, 4,  6, 13, 21, 28);
	/* (a1 + a6 + a9 + a11 + a12 + a15)(b1 + b6 + b9 + b11 + b12 + b15)*(X^4 + X^5 + X^6 + X1 + X2 + X3 + X9 + X^20 + X^21 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 6, a1, a6, a9, a11, a12, a15),
				mzd_sum(t2, 6, b1, b6, b9, b11, b12, b15), 0), X, 12,  4, 5, 6, 11, 12, 13, 19, 20, 21, 26, 27, 28);
	/* (a3 + a6 + a8 + a9 + a12 + a13 + a14 + a15)(b3 + b6 + b8 + b9 + b12 + b13 + b14 + b15)*(X^2 + X^3 + X^5 + X^6 + X^9 + X0 + X2 + X3 + X7 + X8 + X^20 + X^21 + X^24 + X^25 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a3, a6, a8, a9, a12, a13, a14, a15),
				mzd_sum(t2, 8, b3, b6, b8, b9, b12, b13, b14, b15), 0), X, 16,  2, 3, 5, 6, 9, 10, 12, 13, 17, 18, 20, 21, 24, 25, 27, 28);
	/* (a4 + a7 + a9 + a10 + a13 + a14 + a15)(b4 + b7 + b9 + b10 + b13 + b14 + b15)*(X^3 + X^5 + X0 + X2 + X8 + X^20 + X^25 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a4, a7, a9, a10, a13, a14, a15),
				mzd_sum(t2, 7, b4, b7, b9, b10, b13, b14, b15), 0), X, 8,  3, 5, 10, 12, 18, 20, 25, 27);
	/* (a0 + a4 + a5 + a7 + a8 + a9 + a11 + a13)(b0 + b4 + b5 + b7 + b8 + b9 + b11 + b13)*(X^3 + X0 + X8 + X^25) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a4, a5, a7, a8, a9, a11, a13),
				mzd_sum(t2, 8, b0, b4, b5, b7, b8, b9, b11, b13), 0), X, 4,  3, 10, 18, 25);
	/* (a0 + a1 + a5 + a6 + a8 + a9 + a10 + a12 + a14)(b0 + b1 + b5 + b6 + b8 + b9 + b10 + b12 + b14)*(X^2 + X^9 + X7 + X^24) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 9, a0, a1, a5, a6, a8, a9, a10, a12, a14),
				mzd_sum(t2, 9, b0, b1, b5, b6, b8, b9, b10, b12, b14), 0), X, 4,  2, 9, 17, 24);
	/* (a3 + a4 + a6 + a7 + a8 + a10 + a12)(b3 + b4 + b6 + b7 + b8 + b10 + b12)*(X^4 + X1 + X9 + X^26) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a3, a4, a6, a7, a8, a10, a12),
				mzd_sum(t2, 7, b3, b4, b6, b7, b8, b10, b12), 0), X, 4,  4, 11, 19, 26);
	/* (a1 + a2 + a4 + a5 + a6 + a8 + a10 + a15)(b1 + b2 + b4 + b5 + b6 + b8 + b10 + b15)*(X^5 + X^6 + X2 + X3 + X^20 + X^21 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a1, a2, a4, a5, a6, a8, a10, a15),
				mzd_sum(t2, 8, b1, b2, b4, b5, b6, b8, b10, b15), 0), X, 8,  5, 6, 12, 13, 20, 21, 27, 28);
	/* (a0 + a2 + a3 + a6 + a7 + a8 + a9 + a10 + a14 + a15)(b0 + b2 + b3 + b6 + b7 + b8 + b9 + b10 + b14 + b15)*(X^2 + X^3 + X^4 + X^5 + X^9 + X0 + X1 + X2 + X7 + X8 + X9 + X^20 + X^24 + X^25 + X^26 + X^27) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 10, a0, a2, a3, a6, a7, a8, a9, a10, a14, a15),
				mzd_sum(t2, 10, b0, b2, b3, b6, b7, b8, b9, b10, b14, b15), 0), X, 16,  2, 3, 4, 5, 9, 10, 11, 12, 17, 18, 19, 20, 24, 25, 26, 27);
	/* (a0 + a1 + a3 + a4 + a5 + a7 + a9 + a14)(b0 + b1 + b3 + b4 + b5 + b7 + b9 + b14)*(X^2 + X^4 + X^6 + X^9 + X1 + X3 + X7 + X9 + X^21 + X^24 + X^26 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a3, a4, a5, a7, a9, a14),
				mzd_sum(t2, 8, b0, b1, b3, b4, b5, b7, b9, b14), 0), X, 12,  2, 4, 6, 9, 11, 13, 17, 19, 21, 24, 26, 28);
	/* (a0 + a1 + a2 + a4 + a6 + a11 + a14)(b0 + b1 + b2 + b4 + b6 + b11 + b14)*(X^2 + X^4 + X^5 + X^6 + X^9 + X1 + X2 + X3 + X7 + X9 + X^20 + X^21 + X^24 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a2, a4, a6, a11, a14),
				mzd_sum(t2, 7, b0, b1, b2, b4, b6, b11, b14), 0), X, 16,  2, 4, 5, 6, 9, 11, 12, 13, 17, 19, 20, 21, 24, 26, 27, 28);
	/* (a0 + a2 + a3 + a4 + a6 + a8 + a13)(b0 + b2 + b3 + b4 + b6 + b8 + b13)*(X^2 + X^3 + X^4 + X^5 + X^6 + X^9 + X0 + X1 + X2 + X3 + X7 + X8 + X9 + X^20 + X^21 + X^24 + X^25 + X^26 + X^27 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 7, a0, a2, a3, a4, a6, a8, a13),
				mzd_sum(t2, 7, b0, b2, b3, b4, b6, b8, b13), 0), X, 20,  2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 17, 18, 19, 20, 21, 24, 25, 26, 27, 28);
	/* (a0 + a1 + a2 + a3 + a4 + a8 + a9 + a11 + a12 + a13 + a15)(b0 + b1 + b2 + b3 + b4 + b8 + b9 + b11 + b12 + b13 + b15)*(X^3 + X^4 + X^6 + X0 + X1 + X3 + X8 + X9 + X^21 + X^25 + X^26 + X^28) */
	_mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 11, a0, a1, a2, a3, a4, a8, a9, a11, a12, a13, a15),
				mzd_sum(t2, 11, b0, b1, b2, b3, b4, b8, b9, b11, b12, b13, b15), 0), X, 12,  3, 4, 6, 10, 11, 13, 18, 19, 21, 25, 26, 28);

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}
