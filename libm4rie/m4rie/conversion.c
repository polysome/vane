/******************************************************************************
 *
 *            M4RIE: Linear Algebra over GF(2^e)
 *
 *    Copyright (C) 2010-2013 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include "conversion.h"

static inline word word_slice_64_02(word a) {
	a = (a & xcccccccc) | (a & xcccccccc>> 2)<< 1;
	a = (a & xf0f0f0f0) | (a & xf0f0f0f0>> 4)<< 2;
	a = (a & xff00ff00) | (a & xff00ff00>> 8)<< 4;
	a = (a & xffff0000) | (a & xffff0000>>16)<< 8;
	a = (a & xffffffff) | (a & xffffffff>>32)<<16;
	return a;
}

static inline word word_slice_64_04(word a) {
	a = (a & xf0f0f0f0) | (a & xf0f0f0f0>> 4)<< 3;
	a = (a & xff00ff00) | (a & xff00ff00>> 8)<< 6;
	a = (a & xffff0000) | (a & xffff0000>>16)<<12;
	a = (a & xffffffff) | (a & xffffffff>>32)<<24;
	return a;
}

static inline word word_cling_64_02(word a) {
	a = (a & xffff0000 & x__left32) | (a & (xffff0000>>16) & x__left32)>>16;
	a = (a & xff00ff00) | (a & xff00ff00>> 8)>> 8;
	a = (a & xf0f0f0f0) | (a & xf0f0f0f0>> 4)>> 4;
	a = (a & xcccccccc) | (a & xcccccccc>> 2)>> 2;
	a = (a & xaaaaaaaa) | (a & xaaaaaaaa>> 1)>> 1;
	return a;
}

static inline word word_cling_64_04(word a) {
	a = (a & xff00ff00 & x__left16) | (a & (xff00ff00>> 8) & x__left16)>>24;
	a = (a & xf0f0f0f0) | (a & xf0f0f0f0>> 4)>>12;
	a = (a & xcccccccc) | (a & xcccccccc>> 2)>> 6;
	a = (a & xaaaaaaaa) | (a & xaaaaaaaa>> 1)>> 3;
	return a;
}

mzd_slice_t *mzed_slice(mzd_slice_t *A, const mzed_t *Z) {
	if (A == NULL) {
		A = mzd_slice_init(Z->finite_field, Z->nrows, Z->ncols);
	} else {
		mzd_slice_set_ui(A, 0);
	}

	switch(Z->finite_field->degree) {
		case  2: return _mzed_slice2(A,Z);

		case  3: return _mzed_slice4(A,Z);
		case  4: return _mzed_slice4(A,Z);

		case  5: return _mzed_slice8(A,Z);
		case  6: return _mzed_slice8(A,Z);
		case  7: return _mzed_slice8(A,Z);
		case  8: return _mzed_slice8(A,Z);

		case  9: return _mzed_slice16(A,Z);
		case 10: return _mzed_slice16(A,Z);
		case 11: return _mzed_slice16(A,Z);
		case 12: return _mzed_slice16(A,Z);
		case 13: return _mzed_slice16(A,Z);
		case 14: return _mzed_slice16(A,Z);
		case 15: return _mzed_slice16(A,Z);
		case 16: return _mzed_slice16(A,Z);
		default:
				 m4ri_die("slicing not implemented for this degree");
	}
	return A;
}

mzed_t *mzed_cling(mzed_t *A, const mzd_slice_t *Z) {
	if (A == NULL) {
		A = mzed_init(Z->finite_field, Z->nrows, Z->ncols);
	}
	else {
		mzed_set_ui(A, 0);
	}

	switch(Z->finite_field->degree) {
		case  2: return _mzed_cling2(A,Z);

		case  3: return _mzed_cling4(A,Z);
		case  4: return _mzed_cling4(A,Z);

		case  5: return _mzed_cling8(A,Z);
		case  6: return _mzed_cling8(A,Z);
		case  7: return _mzed_cling8(A,Z);
		case  8: return _mzed_cling8(A,Z);

		case  9: return _mzed_cling16(A,Z);
		case 10: return _mzed_cling16(A,Z);
		case 11: return _mzed_cling16(A,Z);
		case 12: return _mzed_cling16(A,Z);
		case 13: return _mzed_cling16(A,Z);
		case 14: return _mzed_cling16(A,Z);
		case 15: return _mzed_cling16(A,Z);
		case 16: return _mzed_cling16(A,Z);
		default:
				 m4ri_die("clinging not implemented for this degree");
	}
	return A;
}

mzd_slice_t *_mzed_slice2(mzd_slice_t *T, const mzed_t *F) {
	assert(T && (T->depth >= 2));
	size_t j, j2 = 0;

	const word bitmask_end = T->x[0]->high_bitmask;
	register word r0,r1,r2,r3;

	if (mzed_is_zero(F))
		return T;

	for(size_t i=0; i<T->nrows; i++) {
		word *t0 = T->x[0]->rows[i];
		word *t1 = T->x[1]->rows[i];
		const word *f  = F->x->rows[i];

		/* bulk of work */
		for(j=0, j2=0; j+2 < F->x->width; j+=2,j2++) {
			r0 = f[j+0], r1 = f[j+1];
			r2 = word_slice_64_02(r0<<1 & xaaaaaaaa);
			r3 = word_slice_64_02(r1<<1 & xaaaaaaaa);
			t0[j2] = r3 | (r2>>32);

			r2 = word_slice_64_02(r0<<0 & xaaaaaaaa);
			r3 = word_slice_64_02(r1<<0 & xaaaaaaaa);
			t1[j2] = r3 | (r2>>32);
		}

		switch(F->x->width - j) {
			case 2:
				r0 = f[j+0]; r1 = f[j+1];

				r2 = word_slice_64_02(r0<<1 & xaaaaaaaa);
				r3 = word_slice_64_02(r1<<1 & xaaaaaaaa);
				t0[j2] &= ~bitmask_end;
				t0[j2] |= (r3 | (r2>>32)) & bitmask_end;

				r2 = word_slice_64_02(r0<<0 & xaaaaaaaa);
				r3 = word_slice_64_02(r1<<0 & xaaaaaaaa);
				t1[j2] &= ~bitmask_end;
				t1[j2] |= (r3 | (r2>>32)) & bitmask_end;
				break;
			case 1:
				r0 = f[j+0];

				r2 = word_slice_64_02(r0<<1 & xaaaaaaaa);
				t0[j2] &= ~bitmask_end;
				t0[j2] |= (r2>>32) & bitmask_end;

				r2 = word_slice_64_02(r0<<0 & xaaaaaaaa);
				t1[j2] &= ~bitmask_end;
				t1[j2] |= (r2>>32) & bitmask_end;
				break;
			default:
				m4ri_die("impossible");
		}
	}

	return T;
}

mzed_t *_mzed_cling2(mzed_t *T, const mzd_slice_t *F) {
	size_t j,j2 = 0;
	register word tmp;

	const word bitmask_end = T->x->high_bitmask;

	if (mzd_slice_is_zero(F))
		return T;

	for(size_t i=0; i<T->nrows; i++) {
		const word *f0 = F->x[0]->rows[i];
		const word *f1 = F->x[1]->rows[i];
		word *t  = T->x->rows[i];

		for(j=0, j2=0; j+2 < T->x->width; j+=2, j2++) {
			t[j+0] = (word_cling_64_02(f0[j2]<<32)>>1) | (word_cling_64_02(f1[j2]<<32)>>0);
			t[j+1] = (word_cling_64_02(f0[j2]<< 0)>>1) | (word_cling_64_02(f1[j2]<< 0)>>0);
		}
		switch(T->x->width - j) {
			case 2:
				tmp    = (word_cling_64_02(f0[j2]<< 0)>>1) | (word_cling_64_02(f1[j2]<< 0)>>0);
				t[j+0] = (word_cling_64_02(f0[j2]<<32)>>1) | (word_cling_64_02(f1[j2]<<32)>>0);
				t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
				break;
			case 1:
				tmp    = (word_cling_64_02(f0[j2]<<32)>>1) | (word_cling_64_02(f1[j2]<<32)>>0);
				t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
				break;
		}
	}
	return T;
}

mzd_slice_t *_mzed_slice4(mzd_slice_t *T, const mzed_t *F) {
	assert(T && (T->depth == 3 || T->depth == 4));
	size_t j, j2 = 0;
	register word r0,r1,r2,r3 = 0;

	const word bitmask_end = T->x[0]->high_bitmask;

	if (mzed_is_zero(F))
		return T;

	if (T->depth == 3) {
		for(size_t i=0; i<T->nrows; i++) {
			word *t0 = T->x[0]->rows[i];
			word *t1 = T->x[1]->rows[i];
			word *t2 = T->x[2]->rows[i];
			const word const *f  = F->x->rows[i];

			/* bulk of work */
			for(j=0, j2=0; j+4 < F->x->width; j+=4,j2++) {
				t0[j2] = word_slice_64_04(f[j+0]<<3 & x88888888)>>48 | word_slice_64_04(f[j+1]<<3 & x88888888)>>32 \
						 |      word_slice_64_04(f[j+2]<<3 & x88888888)>>16 | word_slice_64_04(f[j+3]<<3 & x88888888)>> 0;
				t1[j2] = word_slice_64_04(f[j+0]<<2 & x88888888)>>48 | word_slice_64_04(f[j+1]<<2 & x88888888)>>32 \
						 |      word_slice_64_04(f[j+2]<<2 & x88888888)>>16 | word_slice_64_04(f[j+3]<<2 & x88888888)>> 0;
				t2[j2] = word_slice_64_04(f[j+0]<<1 & x88888888)>>48 | word_slice_64_04(f[j+1]<<1 & x88888888)>>32 \
						 |      word_slice_64_04(f[j+2]<<1 & x88888888)>>16 | word_slice_64_04(f[j+3]<<1 & x88888888)>> 0;
			}
			r0 = r1 = r2 = 0;
			switch(F->x->width - j) {
				case 4:
					r0 |= word_slice_64_04(f[j+3]<<3 & x88888888)>> 0;
					r1 |= word_slice_64_04(f[j+3]<<2 & x88888888)>> 0;
					r2 |= word_slice_64_04(f[j+3]<<1 & x88888888)>> 0;
				case 3:
					r0 |= word_slice_64_04(f[j+2]<<3 & x88888888)>>16;
					r1 |= word_slice_64_04(f[j+2]<<2 & x88888888)>>16;
					r2 |= word_slice_64_04(f[j+2]<<1 & x88888888)>>16;
				case 2:
					r0 |= word_slice_64_04(f[j+1]<<3 & x88888888)>>32;
					r1 |= word_slice_64_04(f[j+1]<<2 & x88888888)>>32;
					r2 |= word_slice_64_04(f[j+1]<<1 & x88888888)>>32;
				case 1:
					r0 |= word_slice_64_04(f[j+0]<<3 & x88888888)>>48;
					r1 |= word_slice_64_04(f[j+0]<<2 & x88888888)>>48;
					r2 |= word_slice_64_04(f[j+0]<<1 & x88888888)>>48;
					break;
				default:
					m4ri_die("impossible");
			}
			t0[j2] |= r0 & bitmask_end;
			t1[j2] |= r1 & bitmask_end;
			t2[j2] |= r2 & bitmask_end;
		}
	} else {
		for(size_t i=0; i<T->nrows; i++) {
			word *t0 = T->x[0]->rows[i];
			word *t1 = T->x[1]->rows[i];
			word *t2 = T->x[2]->rows[i];
			word *t3 = T->x[3]->rows[i];
			const word const *f  = F->x->rows[i];

			/* bulk of work */
			for(j=0, j2=0; j+4 < F->x->width; j+=4,j2++) {
				t0[j2] = word_slice_64_04(f[j+0]<<3 & x88888888)>>48 | word_slice_64_04(f[j+1]<<3 & x88888888)>>32 \
						 |      word_slice_64_04(f[j+2]<<3 & x88888888)>>16 | word_slice_64_04(f[j+3]<<3 & x88888888)>> 0;
				t1[j2] = word_slice_64_04(f[j+0]<<2 & x88888888)>>48 | word_slice_64_04(f[j+1]<<2 & x88888888)>>32 \
						 |      word_slice_64_04(f[j+2]<<2 & x88888888)>>16 | word_slice_64_04(f[j+3]<<2 & x88888888)>> 0;
				t2[j2] = word_slice_64_04(f[j+0]<<1 & x88888888)>>48 | word_slice_64_04(f[j+1]<<1 & x88888888)>>32 \
						 |      word_slice_64_04(f[j+2]<<1 & x88888888)>>16 | word_slice_64_04(f[j+3]<<1 & x88888888)>> 0;
				t3[j2] = word_slice_64_04(f[j+0]<<0 & x88888888)>>48 | word_slice_64_04(f[j+1]<<0 & x88888888)>>32 \
						 |      word_slice_64_04(f[j+2]<<0 & x88888888)>>16 | word_slice_64_04(f[j+3]<<0 & x88888888)>> 0;
			}
			r0 = r1 = r2 = r3 = 0;
			switch(F->x->width - j) {
				case 4:
					r0 |= word_slice_64_04(f[j+3]<<3 & x88888888)>> 0;
					r1 |= word_slice_64_04(f[j+3]<<2 & x88888888)>> 0;
					r2 |= word_slice_64_04(f[j+3]<<1 & x88888888)>> 0;
					r3 |= word_slice_64_04(f[j+3]<<0 & x88888888)>> 0;
				case 3:
					r0 |= word_slice_64_04(f[j+2]<<3 & x88888888)>>16;
					r1 |= word_slice_64_04(f[j+2]<<2 & x88888888)>>16;
					r2 |= word_slice_64_04(f[j+2]<<1 & x88888888)>>16;
					r3 |= word_slice_64_04(f[j+2]<<0 & x88888888)>>16;
				case 2:
					r0 |= word_slice_64_04(f[j+1]<<3 & x88888888)>>32;
					r1 |= word_slice_64_04(f[j+1]<<2 & x88888888)>>32;
					r2 |= word_slice_64_04(f[j+1]<<1 & x88888888)>>32;
					r3 |= word_slice_64_04(f[j+1]<<0 & x88888888)>>32;
				case 1:
					r0 |= word_slice_64_04(f[j+0]<<3 & x88888888)>>48;
					r1 |= word_slice_64_04(f[j+0]<<2 & x88888888)>>48;
					r2 |= word_slice_64_04(f[j+0]<<1 & x88888888)>>48;
					r3 |= word_slice_64_04(f[j+0]<<0 & x88888888)>>48;
					break;
				default:
					m4ri_die("impossible");
			}
			t0[j2] |= r0 & bitmask_end;
			t1[j2] |= r1 & bitmask_end;
			t2[j2] |= r2 & bitmask_end;
			t3[j2] |= r3 & bitmask_end;
		}
	}
	return T;
}

mzed_t *_mzed_cling4(mzed_t *T, const mzd_slice_t *F) {
	size_t j,j2 = 0;

	const word bitmask_end = T->x->high_bitmask;

	if (mzd_slice_is_zero(F))
		return T;

	if (F->finite_field->degree == 4) {
		for(rci_t i=0; i<T->nrows; i++) {
			const word *f0 = F->x[0]->rows[i];
			const word *f1 = F->x[1]->rows[i];
			const word *f2 = F->x[2]->rows[i];
			const word *f3 = F->x[3]->rows[i];
			word *t  = T->x->rows[i];

			for(j=0, j2=0; j+4 < T->x->width; j+=4, j2++) {
				t[j+0] = (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1) | (word_cling_64_04(f3[j2]<<48)>>0);
				t[j+1] = (word_cling_64_04(f0[j2]<<32)>>3) | (word_cling_64_04(f1[j2]<<32)>>2) | (word_cling_64_04(f2[j2]<<32)>>1) | (word_cling_64_04(f3[j2]<<32)>>0);
				t[j+2] = (word_cling_64_04(f0[j2]<<16)>>3) | (word_cling_64_04(f1[j2]<<16)>>2) | (word_cling_64_04(f2[j2]<<16)>>1) | (word_cling_64_04(f3[j2]<<16)>>0);
				t[j+3] = (word_cling_64_04(f0[j2]<< 0)>>3) | (word_cling_64_04(f1[j2]<< 0)>>2) | (word_cling_64_04(f2[j2]<< 0)>>1) | (word_cling_64_04(f3[j2]<< 0)>>0);
			}

			register word tmp=0;
			switch(T->x->width - j) {
				case 4:
					t[j+0] = (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1) | (word_cling_64_04(f3[j2]<<48)>>0);
					t[j+1] = (word_cling_64_04(f0[j2]<<32)>>3) | (word_cling_64_04(f1[j2]<<32)>>2) | (word_cling_64_04(f2[j2]<<32)>>1) | (word_cling_64_04(f3[j2]<<32)>>0);
					t[j+2] = (word_cling_64_04(f0[j2]<<16)>>3) | (word_cling_64_04(f1[j2]<<16)>>2) | (word_cling_64_04(f2[j2]<<16)>>1) | (word_cling_64_04(f3[j2]<<16)>>0);
					tmp    = (word_cling_64_04(f0[j2]<< 0)>>3) | (word_cling_64_04(f1[j2]<< 0)>>2) | (word_cling_64_04(f2[j2]<< 0)>>1) | (word_cling_64_04(f3[j2]<< 0)>>0);
					t[j+3] = (t[j+3] & ~bitmask_end) | (tmp & bitmask_end);
					break;
				case 3:
					t[j+0] = (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1) | (word_cling_64_04(f3[j2]<<48)>>0);
					t[j+1] = (word_cling_64_04(f0[j2]<<32)>>3) | (word_cling_64_04(f1[j2]<<32)>>2) | (word_cling_64_04(f2[j2]<<32)>>1) | (word_cling_64_04(f3[j2]<<32)>>0);
					tmp    = (word_cling_64_04(f0[j2]<<16)>>3) | (word_cling_64_04(f1[j2]<<16)>>2) | (word_cling_64_04(f2[j2]<<16)>>1) | (word_cling_64_04(f3[j2]<<16)>>0);
					t[j+2] = (t[j+2] & ~bitmask_end) | (tmp & bitmask_end);
					break;
				case 2:
					t[j+0] = (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1) | (word_cling_64_04(f3[j2]<<48)>>0);
					tmp    = (word_cling_64_04(f0[j2]<<32)>>3) | (word_cling_64_04(f1[j2]<<32)>>2) | (word_cling_64_04(f2[j2]<<32)>>1) | (word_cling_64_04(f3[j2]<<32)>>0);
					t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
					break;
				case 1:
					tmp =    (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1) | (word_cling_64_04(f3[j2]<<48)>>0);
					t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
					break;
				default:
					m4ri_die("impossible");
			}
		}
	} else { //degree == 3
		for(rci_t i=0; i<T->nrows; i++) {
			const word *f0 = F->x[0]->rows[i];
			const word *f1 = F->x[1]->rows[i];
			const word *f2 = F->x[2]->rows[i];
			word *t  = T->x->rows[i];

			for(j=0, j2=0; j+4 < T->x->width; j+=4, j2++) {
				t[j+0] = (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1);
				t[j+1] = (word_cling_64_04(f0[j2]<<32)>>3) | (word_cling_64_04(f1[j2]<<32)>>2) | (word_cling_64_04(f2[j2]<<32)>>1);
				t[j+2] = (word_cling_64_04(f0[j2]<<16)>>3) | (word_cling_64_04(f1[j2]<<16)>>2) | (word_cling_64_04(f2[j2]<<16)>>1);
				t[j+3] = (word_cling_64_04(f0[j2]<< 0)>>3) | (word_cling_64_04(f1[j2]<< 0)>>2) | (word_cling_64_04(f2[j2]<< 0)>>1);
			}

			register word tmp=0;
			switch(T->x->width - j) {
				case 4:
					t[j+0] = (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1);
					t[j+1] = (word_cling_64_04(f0[j2]<<32)>>3) | (word_cling_64_04(f1[j2]<<32)>>2) | (word_cling_64_04(f2[j2]<<32)>>1);
					t[j+2] = (word_cling_64_04(f0[j2]<<16)>>3) | (word_cling_64_04(f1[j2]<<16)>>2) | (word_cling_64_04(f2[j2]<<16)>>1);
					tmp  =   (word_cling_64_04(f0[j2]<< 0)>>3) | (word_cling_64_04(f1[j2]<< 0)>>2) | (word_cling_64_04(f2[j2]<< 0)>>1);
					t[j+3] = (t[j+3] & ~bitmask_end) | (tmp & bitmask_end);
					break;
				case 3:
					t[j+0] = (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1);
					t[j+1] = (word_cling_64_04(f0[j2]<<32)>>3) | (word_cling_64_04(f1[j2]<<32)>>2) | (word_cling_64_04(f2[j2]<<32)>>1);
					tmp =    (word_cling_64_04(f0[j2]<<16)>>3) | (word_cling_64_04(f1[j2]<<16)>>2) | (word_cling_64_04(f2[j2]<<16)>>1);
					t[j+2] = (t[j+2] & ~bitmask_end) | (tmp & bitmask_end);
					break;
				case 2:
					t[j+0] = (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1);
					tmp =    (word_cling_64_04(f0[j2]<<32)>>3) | (word_cling_64_04(f1[j2]<<32)>>2) | (word_cling_64_04(f2[j2]<<32)>>1);
					t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
					break;
				case 1:
					tmp =    (word_cling_64_04(f0[j2]<<48)>>3) | (word_cling_64_04(f1[j2]<<48)>>2) | (word_cling_64_04(f2[j2]<<48)>>1);
					t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
					break;
				default:
					m4ri_die("impossible");
			}
		}
	}
	return T;
}


