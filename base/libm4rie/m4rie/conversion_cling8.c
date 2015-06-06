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

static inline word word_cling_64_08(word a) {
	a = (a & xf0f0f0f0 & x__left08) | (a & xf0f0f0f0>> 4  & x__left08)>>28;
	a = (a & xcccccccc) | (a & xcccccccc>> 2)>>14;
	a = (a & xaaaaaaaa) | (a & xaaaaaaaa>> 1)>> 7;
	return a;
}

mzed_t *_mzed_cling8(mzed_t *T, const mzd_slice_t *F) {
	size_t j,j2 = 0;

	const word bitmask_end = T->x->high_bitmask;

	if (mzd_slice_is_zero(F))
		return T;

	switch (F->finite_field->degree) {
		case 8: {
					for(rci_t i=0; i<T->nrows; i++) {
						const word *f0 = F->x[0]->rows[i];
						const word *f1 = F->x[1]->rows[i];
						const word *f2 = F->x[2]->rows[i];
						const word *f3 = F->x[3]->rows[i];
						const word *f4 = F->x[4]->rows[i];
						const word *f5 = F->x[5]->rows[i];
						const word *f6 = F->x[6]->rows[i];
						const word *f7 = F->x[7]->rows[i];
						word *t  = T->x->rows[i];

						for(j=0, j2=0; j+8 < T->x->width; j+=8, j2++) {
							t[j+0] = (word_cling_64_08(f0[j2]<<56)>>7) | (word_cling_64_08(f1[j2]<<56)>>6) | (word_cling_64_08(f2[j2]<<56)>>5) | (word_cling_64_08(f3[j2]<<56)>>4) \
									 |      (word_cling_64_08(f4[j2]<<56)>>3) | (word_cling_64_08(f5[j2]<<56)>>2) | (word_cling_64_08(f6[j2]<<56)>>1) | (word_cling_64_08(f7[j2]<<56)>>0);
							t[j+1] = (word_cling_64_08(f0[j2]<<48)>>7) | (word_cling_64_08(f1[j2]<<48)>>6) | (word_cling_64_08(f2[j2]<<48)>>5) | (word_cling_64_08(f3[j2]<<48)>>4) \
									 |      (word_cling_64_08(f4[j2]<<48)>>3) | (word_cling_64_08(f5[j2]<<48)>>2) | (word_cling_64_08(f6[j2]<<48)>>1) | (word_cling_64_08(f7[j2]<<48)>>0);
							t[j+2] = (word_cling_64_08(f0[j2]<<40)>>7) | (word_cling_64_08(f1[j2]<<40)>>6) | (word_cling_64_08(f2[j2]<<40)>>5) | (word_cling_64_08(f3[j2]<<40)>>4) \
									 |      (word_cling_64_08(f4[j2]<<40)>>3) | (word_cling_64_08(f5[j2]<<40)>>2) | (word_cling_64_08(f6[j2]<<40)>>1) | (word_cling_64_08(f7[j2]<<40)>>0);
							t[j+3] = (word_cling_64_08(f0[j2]<<32)>>7) | (word_cling_64_08(f1[j2]<<32)>>6) | (word_cling_64_08(f2[j2]<<32)>>5) | (word_cling_64_08(f3[j2]<<32)>>4) \
									 |      (word_cling_64_08(f4[j2]<<32)>>3) | (word_cling_64_08(f5[j2]<<32)>>2) | (word_cling_64_08(f6[j2]<<32)>>1) | (word_cling_64_08(f7[j2]<<32)>>0);
							t[j+4] = (word_cling_64_08(f0[j2]<<24)>>7) | (word_cling_64_08(f1[j2]<<24)>>6) | (word_cling_64_08(f2[j2]<<24)>>5) | (word_cling_64_08(f3[j2]<<24)>>4) \
									 |      (word_cling_64_08(f4[j2]<<24)>>3) | (word_cling_64_08(f5[j2]<<24)>>2) | (word_cling_64_08(f6[j2]<<24)>>1) | (word_cling_64_08(f7[j2]<<24)>>0);
							t[j+5] = (word_cling_64_08(f0[j2]<<16)>>7) | (word_cling_64_08(f1[j2]<<16)>>6) | (word_cling_64_08(f2[j2]<<16)>>5) | (word_cling_64_08(f3[j2]<<16)>>4) \
									 |      (word_cling_64_08(f4[j2]<<16)>>3) | (word_cling_64_08(f5[j2]<<16)>>2) | (word_cling_64_08(f6[j2]<<16)>>1) | (word_cling_64_08(f7[j2]<<16)>>0);
							t[j+6] = (word_cling_64_08(f0[j2]<< 8)>>7) | (word_cling_64_08(f1[j2]<< 8)>>6) | (word_cling_64_08(f2[j2]<< 8)>>5) | (word_cling_64_08(f3[j2]<< 8)>>4) \
									 |      (word_cling_64_08(f4[j2]<< 8)>>3) | (word_cling_64_08(f5[j2]<< 8)>>2) | (word_cling_64_08(f6[j2]<< 8)>>1) | (word_cling_64_08(f7[j2]<< 8)>>0);
							t[j+7] = (word_cling_64_08(f0[j2]<< 0)>>7) | (word_cling_64_08(f1[j2]<< 0)>>6) | (word_cling_64_08(f2[j2]<< 0)>>5) | (word_cling_64_08(f3[j2]<< 0)>>4) \
									 |      (word_cling_64_08(f4[j2]<< 0)>>3) | (word_cling_64_08(f5[j2]<< 0)>>2) | (word_cling_64_08(f6[j2]<< 0)>>1) | (word_cling_64_08(f7[j2]<< 0)>>0);
						}

						register word tmp = t[T->x->width-1];
						switch(T->x->width - j) {
							case 8:
								t[j+7] = (word_cling_64_08(f0[j2]<< 0)>>7) | (word_cling_64_08(f1[j2]<< 0)>>6) | (word_cling_64_08(f2[j2]<< 0)>>5) | (word_cling_64_08(f3[j2]<< 0)>>4) \
										 |      (word_cling_64_08(f4[j2]<< 0)>>3) | (word_cling_64_08(f5[j2]<< 0)>>2) | (word_cling_64_08(f6[j2]<< 0)>>1) | (word_cling_64_08(f7[j2]<< 0)>>0);
							case 7:
								t[j+6] = (word_cling_64_08(f0[j2]<< 8)>>7) | (word_cling_64_08(f1[j2]<< 8)>>6) | (word_cling_64_08(f2[j2]<< 8)>>5) | (word_cling_64_08(f3[j2]<< 8)>>4) \
										 |      (word_cling_64_08(f4[j2]<< 8)>>3) | (word_cling_64_08(f5[j2]<< 8)>>2) | (word_cling_64_08(f6[j2]<< 8)>>1) | (word_cling_64_08(f7[j2]<< 8)>>0);
							case 6:
								t[j+5] = (word_cling_64_08(f0[j2]<<16)>>7) | (word_cling_64_08(f1[j2]<<16)>>6) | (word_cling_64_08(f2[j2]<<16)>>5) | (word_cling_64_08(f3[j2]<<16)>>4) \
										 |      (word_cling_64_08(f4[j2]<<16)>>3) | (word_cling_64_08(f5[j2]<<16)>>2) | (word_cling_64_08(f6[j2]<<16)>>1) | (word_cling_64_08(f7[j2]<<16)>>0);
							case 5:
								t[j+4] = (word_cling_64_08(f0[j2]<<24)>>7) | (word_cling_64_08(f1[j2]<<24)>>6) | (word_cling_64_08(f2[j2]<<24)>>5) | (word_cling_64_08(f3[j2]<<24)>>4) \
										 |      (word_cling_64_08(f4[j2]<<24)>>3) | (word_cling_64_08(f5[j2]<<24)>>2) | (word_cling_64_08(f6[j2]<<24)>>1) | (word_cling_64_08(f7[j2]<<24)>>0);
							case 4:
								t[j+3] = (word_cling_64_08(f0[j2]<<32)>>7) | (word_cling_64_08(f1[j2]<<32)>>6) | (word_cling_64_08(f2[j2]<<32)>>5) | (word_cling_64_08(f3[j2]<<32)>>4) \
										 |      (word_cling_64_08(f4[j2]<<32)>>3) | (word_cling_64_08(f5[j2]<<32)>>2) | (word_cling_64_08(f6[j2]<<32)>>1) | (word_cling_64_08(f7[j2]<<32)>>0);
							case 3:
								t[j+2] = (word_cling_64_08(f0[j2]<<40)>>7) | (word_cling_64_08(f1[j2]<<40)>>6) | (word_cling_64_08(f2[j2]<<40)>>5) | (word_cling_64_08(f3[j2]<<40)>>4) \
										 |      (word_cling_64_08(f4[j2]<<40)>>3) | (word_cling_64_08(f5[j2]<<40)>>2) | (word_cling_64_08(f6[j2]<<40)>>1) | (word_cling_64_08(f7[j2]<<40)>>0);
							case 2:
								t[j+1] = (word_cling_64_08(f0[j2]<<48)>>7) | (word_cling_64_08(f1[j2]<<48)>>6) | (word_cling_64_08(f2[j2]<<48)>>5) | (word_cling_64_08(f3[j2]<<48)>>4) \
										 |      (word_cling_64_08(f4[j2]<<48)>>3) | (word_cling_64_08(f5[j2]<<48)>>2) | (word_cling_64_08(f6[j2]<<48)>>1) | (word_cling_64_08(f7[j2]<<48)>>0);
							case 1:
								t[j+0] = (word_cling_64_08(f0[j2]<<56)>>7) | (word_cling_64_08(f1[j2]<<56)>>6) | (word_cling_64_08(f2[j2]<<56)>>5) | (word_cling_64_08(f3[j2]<<56)>>4) \
										 |      (word_cling_64_08(f4[j2]<<56)>>3) | (word_cling_64_08(f5[j2]<<56)>>2) | (word_cling_64_08(f6[j2]<<56)>>1) | (word_cling_64_08(f7[j2]<<56)>>0);
								break;
							default:
								m4ri_die("impossible");
						}
						t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
					} // for loop
				}
				break;

		case 7: {
					for(rci_t i=0; i<T->nrows; i++) {
						const word *f0 = F->x[0]->rows[i];
						const word *f1 = F->x[1]->rows[i];
						const word *f2 = F->x[2]->rows[i];
						const word *f3 = F->x[3]->rows[i];
						const word *f4 = F->x[4]->rows[i];
						const word *f5 = F->x[5]->rows[i];
						const word *f6 = F->x[6]->rows[i];
						word *t  = T->x->rows[i];

						for(j=0, j2=0; j+8 < T->x->width; j+=8, j2++) {
							t[j+0] = (word_cling_64_08(f0[j2]<<56)>>7) | (word_cling_64_08(f1[j2]<<56)>>6) | (word_cling_64_08(f2[j2]<<56)>>5) | (word_cling_64_08(f3[j2]<<56)>>4) \
									 |      (word_cling_64_08(f4[j2]<<56)>>3) | (word_cling_64_08(f5[j2]<<56)>>2) | (word_cling_64_08(f6[j2]<<56)>>1);
							t[j+1] = (word_cling_64_08(f0[j2]<<48)>>7) | (word_cling_64_08(f1[j2]<<48)>>6) | (word_cling_64_08(f2[j2]<<48)>>5) | (word_cling_64_08(f3[j2]<<48)>>4) \
									 |      (word_cling_64_08(f4[j2]<<48)>>3) | (word_cling_64_08(f5[j2]<<48)>>2) | (word_cling_64_08(f6[j2]<<48)>>1);
							t[j+2] = (word_cling_64_08(f0[j2]<<40)>>7) | (word_cling_64_08(f1[j2]<<40)>>6) | (word_cling_64_08(f2[j2]<<40)>>5) | (word_cling_64_08(f3[j2]<<40)>>4) \
									 |      (word_cling_64_08(f4[j2]<<40)>>3) | (word_cling_64_08(f5[j2]<<40)>>2) | (word_cling_64_08(f6[j2]<<40)>>1);
							t[j+3] = (word_cling_64_08(f0[j2]<<32)>>7) | (word_cling_64_08(f1[j2]<<32)>>6) | (word_cling_64_08(f2[j2]<<32)>>5) | (word_cling_64_08(f3[j2]<<32)>>4) \
									 |      (word_cling_64_08(f4[j2]<<32)>>3) | (word_cling_64_08(f5[j2]<<32)>>2) | (word_cling_64_08(f6[j2]<<32)>>1);
							t[j+4] = (word_cling_64_08(f0[j2]<<24)>>7) | (word_cling_64_08(f1[j2]<<24)>>6) | (word_cling_64_08(f2[j2]<<24)>>5) | (word_cling_64_08(f3[j2]<<24)>>4) \
									 |      (word_cling_64_08(f4[j2]<<24)>>3) | (word_cling_64_08(f5[j2]<<24)>>2) | (word_cling_64_08(f6[j2]<<24)>>1);
							t[j+5] = (word_cling_64_08(f0[j2]<<16)>>7) | (word_cling_64_08(f1[j2]<<16)>>6) | (word_cling_64_08(f2[j2]<<16)>>5) | (word_cling_64_08(f3[j2]<<16)>>4) \
									 |      (word_cling_64_08(f4[j2]<<16)>>3) | (word_cling_64_08(f5[j2]<<16)>>2) | (word_cling_64_08(f6[j2]<<16)>>1);
							t[j+6] = (word_cling_64_08(f0[j2]<< 8)>>7) | (word_cling_64_08(f1[j2]<< 8)>>6) | (word_cling_64_08(f2[j2]<< 8)>>5) | (word_cling_64_08(f3[j2]<< 8)>>4) \
									 |      (word_cling_64_08(f4[j2]<< 8)>>3) | (word_cling_64_08(f5[j2]<< 8)>>2) | (word_cling_64_08(f6[j2]<< 8)>>1);
							t[j+7] = (word_cling_64_08(f0[j2]<< 0)>>7) | (word_cling_64_08(f1[j2]<< 0)>>6) | (word_cling_64_08(f2[j2]<< 0)>>5) | (word_cling_64_08(f3[j2]<< 0)>>4) \
									 |      (word_cling_64_08(f4[j2]<< 0)>>3) | (word_cling_64_08(f5[j2]<< 0)>>2) | (word_cling_64_08(f6[j2]<< 0)>>1);
						}

						register word tmp= t[T->x->width-1];
						switch(T->x->width - j) {
							case 8:
								t[j+7] = (word_cling_64_08(f0[j2]<< 0)>>7) | (word_cling_64_08(f1[j2]<< 0)>>6) | (word_cling_64_08(f2[j2]<< 0)>>5) | (word_cling_64_08(f3[j2]<< 0)>>4) \
										 |      (word_cling_64_08(f4[j2]<< 0)>>3) | (word_cling_64_08(f5[j2]<< 0)>>2) | (word_cling_64_08(f6[j2]<< 0)>>1);
							case 7:
								t[j+6] = (word_cling_64_08(f0[j2]<< 8)>>7) | (word_cling_64_08(f1[j2]<< 8)>>6) | (word_cling_64_08(f2[j2]<< 8)>>5) | (word_cling_64_08(f3[j2]<< 8)>>4) \
										 |      (word_cling_64_08(f4[j2]<< 8)>>3) | (word_cling_64_08(f5[j2]<< 8)>>2) | (word_cling_64_08(f6[j2]<< 8)>>1);
							case 6:
								t[j+5] = (word_cling_64_08(f0[j2]<<16)>>7) | (word_cling_64_08(f1[j2]<<16)>>6) | (word_cling_64_08(f2[j2]<<16)>>5) | (word_cling_64_08(f3[j2]<<16)>>4) \
										 |      (word_cling_64_08(f4[j2]<<16)>>3) | (word_cling_64_08(f5[j2]<<16)>>2) | (word_cling_64_08(f6[j2]<<16)>>1);
							case 5:
								t[j+4] = (word_cling_64_08(f0[j2]<<24)>>7) | (word_cling_64_08(f1[j2]<<24)>>6) | (word_cling_64_08(f2[j2]<<24)>>5) | (word_cling_64_08(f3[j2]<<24)>>4) \
										 |      (word_cling_64_08(f4[j2]<<24)>>3) | (word_cling_64_08(f5[j2]<<24)>>2) | (word_cling_64_08(f6[j2]<<24)>>1);
							case 4:
								t[j+3] = (word_cling_64_08(f0[j2]<<32)>>7) | (word_cling_64_08(f1[j2]<<32)>>6) | (word_cling_64_08(f2[j2]<<32)>>5) | (word_cling_64_08(f3[j2]<<32)>>4) \
										 |      (word_cling_64_08(f4[j2]<<32)>>3) | (word_cling_64_08(f5[j2]<<32)>>2) | (word_cling_64_08(f6[j2]<<32)>>1);
							case 3:
								t[j+2] = (word_cling_64_08(f0[j2]<<40)>>7) | (word_cling_64_08(f1[j2]<<40)>>6) | (word_cling_64_08(f2[j2]<<40)>>5) | (word_cling_64_08(f3[j2]<<40)>>4) \
										 |      (word_cling_64_08(f4[j2]<<40)>>3) | (word_cling_64_08(f5[j2]<<40)>>2) | (word_cling_64_08(f6[j2]<<40)>>1);
							case 2:
								t[j+1] = (word_cling_64_08(f0[j2]<<48)>>7) | (word_cling_64_08(f1[j2]<<48)>>6) | (word_cling_64_08(f2[j2]<<48)>>5) | (word_cling_64_08(f3[j2]<<48)>>4) \
										 |      (word_cling_64_08(f4[j2]<<48)>>3) | (word_cling_64_08(f5[j2]<<48)>>2) | (word_cling_64_08(f6[j2]<<48)>>1);
							case 1:
								t[j+0] = (word_cling_64_08(f0[j2]<<56)>>7) | (word_cling_64_08(f1[j2]<<56)>>6) | (word_cling_64_08(f2[j2]<<56)>>5) | (word_cling_64_08(f3[j2]<<56)>>4) \
										 |      (word_cling_64_08(f4[j2]<<56)>>3) | (word_cling_64_08(f5[j2]<<56)>>2) | (word_cling_64_08(f6[j2]<<56)>>1);
								break;
							default:
								m4ri_die("impossible");
						}
						t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
					}
				}
				break;

		case 6: {
					for(rci_t i=0; i<T->nrows; i++) {
						const word *f0 = F->x[0]->rows[i];
						const word *f1 = F->x[1]->rows[i];
						const word *f2 = F->x[2]->rows[i];
						const word *f3 = F->x[3]->rows[i];
						const word *f4 = F->x[4]->rows[i];
						const word *f5 = F->x[5]->rows[i];
						word *t  = T->x->rows[i];

						for(j=0, j2=0; j+8 < T->x->width; j+=8, j2++) {
							t[j+0] = (word_cling_64_08(f0[j2]<<56)>>7) | (word_cling_64_08(f1[j2]<<56)>>6) | (word_cling_64_08(f2[j2]<<56)>>5) | (word_cling_64_08(f3[j2]<<56)>>4) \
									 |      (word_cling_64_08(f4[j2]<<56)>>3) | (word_cling_64_08(f5[j2]<<56)>>2);
							t[j+1] = (word_cling_64_08(f0[j2]<<48)>>7) | (word_cling_64_08(f1[j2]<<48)>>6) | (word_cling_64_08(f2[j2]<<48)>>5) | (word_cling_64_08(f3[j2]<<48)>>4) \
									 |      (word_cling_64_08(f4[j2]<<48)>>3) | (word_cling_64_08(f5[j2]<<48)>>2);
							t[j+2] = (word_cling_64_08(f0[j2]<<40)>>7) | (word_cling_64_08(f1[j2]<<40)>>6) | (word_cling_64_08(f2[j2]<<40)>>5) | (word_cling_64_08(f3[j2]<<40)>>4) \
									 |      (word_cling_64_08(f4[j2]<<40)>>3) | (word_cling_64_08(f5[j2]<<40)>>2);
							t[j+3] = (word_cling_64_08(f0[j2]<<32)>>7) | (word_cling_64_08(f1[j2]<<32)>>6) | (word_cling_64_08(f2[j2]<<32)>>5) | (word_cling_64_08(f3[j2]<<32)>>4) \
									 |      (word_cling_64_08(f4[j2]<<32)>>3) | (word_cling_64_08(f5[j2]<<32)>>2);
							t[j+4] = (word_cling_64_08(f0[j2]<<24)>>7) | (word_cling_64_08(f1[j2]<<24)>>6) | (word_cling_64_08(f2[j2]<<24)>>5) | (word_cling_64_08(f3[j2]<<24)>>4) \
									 |      (word_cling_64_08(f4[j2]<<24)>>3) | (word_cling_64_08(f5[j2]<<24)>>2);
							t[j+5] = (word_cling_64_08(f0[j2]<<16)>>7) | (word_cling_64_08(f1[j2]<<16)>>6) | (word_cling_64_08(f2[j2]<<16)>>5) | (word_cling_64_08(f3[j2]<<16)>>4) \
									 |      (word_cling_64_08(f4[j2]<<16)>>3) | (word_cling_64_08(f5[j2]<<16)>>2);
							t[j+6] = (word_cling_64_08(f0[j2]<< 8)>>7) | (word_cling_64_08(f1[j2]<< 8)>>6) | (word_cling_64_08(f2[j2]<< 8)>>5) | (word_cling_64_08(f3[j2]<< 8)>>4) \
									 |      (word_cling_64_08(f4[j2]<< 8)>>3) | (word_cling_64_08(f5[j2]<< 8)>>2);
							t[j+7] = (word_cling_64_08(f0[j2]<< 0)>>7) | (word_cling_64_08(f1[j2]<< 0)>>6) | (word_cling_64_08(f2[j2]<< 0)>>5) | (word_cling_64_08(f3[j2]<< 0)>>4) \
									 |      (word_cling_64_08(f4[j2]<< 0)>>3) | (word_cling_64_08(f5[j2]<< 0)>>2);
						}

						register word tmp = t[T->x->width-1];
						switch(T->x->width - j) {
							case 8:
								t[j+7] = (word_cling_64_08(f0[j2]<< 0)>>7) | (word_cling_64_08(f1[j2]<< 0)>>6) | (word_cling_64_08(f2[j2]<< 0)>>5) | (word_cling_64_08(f3[j2]<< 0)>>4) \
										 |      (word_cling_64_08(f4[j2]<< 0)>>3) | (word_cling_64_08(f5[j2]<< 0)>>2);
							case 7:
								t[j+6] = (word_cling_64_08(f0[j2]<< 8)>>7) | (word_cling_64_08(f1[j2]<< 8)>>6) | (word_cling_64_08(f2[j2]<< 8)>>5) | (word_cling_64_08(f3[j2]<< 8)>>4) \
										 |      (word_cling_64_08(f4[j2]<< 8)>>3) | (word_cling_64_08(f5[j2]<< 8)>>2);
							case 6:
								t[j+5] = (word_cling_64_08(f0[j2]<<16)>>7) | (word_cling_64_08(f1[j2]<<16)>>6) | (word_cling_64_08(f2[j2]<<16)>>5) | (word_cling_64_08(f3[j2]<<16)>>4) \
										 |      (word_cling_64_08(f4[j2]<<16)>>3) | (word_cling_64_08(f5[j2]<<16)>>2);
							case 5:
								t[j+4] = (word_cling_64_08(f0[j2]<<24)>>7) | (word_cling_64_08(f1[j2]<<24)>>6) | (word_cling_64_08(f2[j2]<<24)>>5) | (word_cling_64_08(f3[j2]<<24)>>4) \
										 |      (word_cling_64_08(f4[j2]<<24)>>3) | (word_cling_64_08(f5[j2]<<24)>>2);
							case 4:
								t[j+3] = (word_cling_64_08(f0[j2]<<32)>>7) | (word_cling_64_08(f1[j2]<<32)>>6) | (word_cling_64_08(f2[j2]<<32)>>5) | (word_cling_64_08(f3[j2]<<32)>>4) \
										 |      (word_cling_64_08(f4[j2]<<32)>>3) | (word_cling_64_08(f5[j2]<<32)>>2);
							case 3:
								t[j+2] = (word_cling_64_08(f0[j2]<<40)>>7) | (word_cling_64_08(f1[j2]<<40)>>6) | (word_cling_64_08(f2[j2]<<40)>>5) | (word_cling_64_08(f3[j2]<<40)>>4) \
										 |      (word_cling_64_08(f4[j2]<<40)>>3) | (word_cling_64_08(f5[j2]<<40)>>2);
							case 2:
								t[j+1] = (word_cling_64_08(f0[j2]<<48)>>7) | (word_cling_64_08(f1[j2]<<48)>>6) | (word_cling_64_08(f2[j2]<<48)>>5) | (word_cling_64_08(f3[j2]<<48)>>4) \
										 |      (word_cling_64_08(f4[j2]<<48)>>3) | (word_cling_64_08(f5[j2]<<48)>>2);
							case 1:
								t[j+0] = (word_cling_64_08(f0[j2]<<56)>>7) | (word_cling_64_08(f1[j2]<<56)>>6) | (word_cling_64_08(f2[j2]<<56)>>5) | (word_cling_64_08(f3[j2]<<56)>>4) \
										 |      (word_cling_64_08(f4[j2]<<56)>>3) | (word_cling_64_08(f5[j2]<<56)>>2);
								break;
							default:
								m4ri_die("impossible");
						}
						t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
					}
				}
				break;

		case 5: {
					for(rci_t i=0; i<T->nrows; i++) {
						const word *f0 = F->x[0]->rows[i];
						const word *f1 = F->x[1]->rows[i];
						const word *f2 = F->x[2]->rows[i];
						const word *f3 = F->x[3]->rows[i];
						const word *f4 = F->x[4]->rows[i];
						word *t  = T->x->rows[i];

						for(j=0, j2=0; j+8 < T->x->width; j+=8, j2++) {
							t[j+0] = (word_cling_64_08(f0[j2]<<56)>>7) | (word_cling_64_08(f1[j2]<<56)>>6) | (word_cling_64_08(f2[j2]<<56)>>5) | (word_cling_64_08(f3[j2]<<56)>>4) | (word_cling_64_08(f4[j2]<<56)>>3);
							t[j+1] = (word_cling_64_08(f0[j2]<<48)>>7) | (word_cling_64_08(f1[j2]<<48)>>6) | (word_cling_64_08(f2[j2]<<48)>>5) | (word_cling_64_08(f3[j2]<<48)>>4) | (word_cling_64_08(f4[j2]<<48)>>3);
							t[j+2] = (word_cling_64_08(f0[j2]<<40)>>7) | (word_cling_64_08(f1[j2]<<40)>>6) | (word_cling_64_08(f2[j2]<<40)>>5) | (word_cling_64_08(f3[j2]<<40)>>4) | (word_cling_64_08(f4[j2]<<40)>>3);
							t[j+3] = (word_cling_64_08(f0[j2]<<32)>>7) | (word_cling_64_08(f1[j2]<<32)>>6) | (word_cling_64_08(f2[j2]<<32)>>5) | (word_cling_64_08(f3[j2]<<32)>>4) | (word_cling_64_08(f4[j2]<<32)>>3);
							t[j+4] = (word_cling_64_08(f0[j2]<<24)>>7) | (word_cling_64_08(f1[j2]<<24)>>6) | (word_cling_64_08(f2[j2]<<24)>>5) | (word_cling_64_08(f3[j2]<<24)>>4) | (word_cling_64_08(f4[j2]<<24)>>3);
							t[j+5] = (word_cling_64_08(f0[j2]<<16)>>7) | (word_cling_64_08(f1[j2]<<16)>>6) | (word_cling_64_08(f2[j2]<<16)>>5) | (word_cling_64_08(f3[j2]<<16)>>4) | (word_cling_64_08(f4[j2]<<16)>>3);
							t[j+6] = (word_cling_64_08(f0[j2]<< 8)>>7) | (word_cling_64_08(f1[j2]<< 8)>>6) | (word_cling_64_08(f2[j2]<< 8)>>5) | (word_cling_64_08(f3[j2]<< 8)>>4) | (word_cling_64_08(f4[j2]<< 8)>>3);
							t[j+7] = (word_cling_64_08(f0[j2]<< 0)>>7) | (word_cling_64_08(f1[j2]<< 0)>>6) | (word_cling_64_08(f2[j2]<< 0)>>5) | (word_cling_64_08(f3[j2]<< 0)>>4) | (word_cling_64_08(f4[j2]<< 0)>>3);
						}

						register word tmp = t[T->x->width - 1];
						switch(T->x->width - j) {
							case 8: t[j+7] = (word_cling_64_08(f0[j2]<< 0)>>7) | (word_cling_64_08(f1[j2]<< 0)>>6) | (word_cling_64_08(f2[j2]<< 0)>>5) | (word_cling_64_08(f3[j2]<< 0)>>4) | (word_cling_64_08(f4[j2]<< 0)>>3);
							case 7: t[j+6] = (word_cling_64_08(f0[j2]<< 8)>>7) | (word_cling_64_08(f1[j2]<< 8)>>6) | (word_cling_64_08(f2[j2]<< 8)>>5) | (word_cling_64_08(f3[j2]<< 8)>>4) | (word_cling_64_08(f4[j2]<< 8)>>3);
							case 6: t[j+5] = (word_cling_64_08(f0[j2]<<16)>>7) | (word_cling_64_08(f1[j2]<<16)>>6) | (word_cling_64_08(f2[j2]<<16)>>5) | (word_cling_64_08(f3[j2]<<16)>>4) | (word_cling_64_08(f4[j2]<<16)>>3);
							case 5: t[j+4] = (word_cling_64_08(f0[j2]<<24)>>7) | (word_cling_64_08(f1[j2]<<24)>>6) | (word_cling_64_08(f2[j2]<<24)>>5) | (word_cling_64_08(f3[j2]<<24)>>4) | (word_cling_64_08(f4[j2]<<24)>>3);
							case 4: t[j+3] = (word_cling_64_08(f0[j2]<<32)>>7) | (word_cling_64_08(f1[j2]<<32)>>6) | (word_cling_64_08(f2[j2]<<32)>>5) | (word_cling_64_08(f3[j2]<<32)>>4) | (word_cling_64_08(f4[j2]<<32)>>3);
							case 3: t[j+2] = (word_cling_64_08(f0[j2]<<40)>>7) | (word_cling_64_08(f1[j2]<<40)>>6) | (word_cling_64_08(f2[j2]<<40)>>5) | (word_cling_64_08(f3[j2]<<40)>>4) | (word_cling_64_08(f4[j2]<<40)>>3);
							case 2: t[j+1] = (word_cling_64_08(f0[j2]<<48)>>7) | (word_cling_64_08(f1[j2]<<48)>>6) | (word_cling_64_08(f2[j2]<<48)>>5) | (word_cling_64_08(f3[j2]<<48)>>4) | (word_cling_64_08(f4[j2]<<48)>>3);
							case 1: t[j+0] = (word_cling_64_08(f0[j2]<<56)>>7) | (word_cling_64_08(f1[j2]<<56)>>6) | (word_cling_64_08(f2[j2]<<56)>>5) | (word_cling_64_08(f3[j2]<<56)>>4) | (word_cling_64_08(f4[j2]<<56)>>3);
									break;
							default:
									m4ri_die("impossible");
						}
						t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
					}
				}
				break;
		default:
				m4ri_die("impossible");
	}
	return T;
}
