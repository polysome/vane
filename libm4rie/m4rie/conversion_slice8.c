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

static inline word word_slice_64_08(word a) {
	a = (a & xff00ff00) | (a & xff00ff00>> 8)<< 7;
	a = (a & xffff0000) | (a & xffff0000>>16)<<14;
	a = (a & xffffffff) | (a & xffffffff>>32)<<28;
	return a;
}

mzd_slice_t *_mzed_slice8(mzd_slice_t *T, const mzed_t *F) {
	assert(T && (4 < T->depth && T->depth <= 8));
	size_t j, j2 = 0;
	register word r0,r1,r2,r3,r4,r5,r6,r7 = 0;

	const word bitmask_end = T->x[0]->high_bitmask;

	if (mzed_is_zero(F))
		return T;

	switch(T->depth) {
		case 8: {
					for(size_t i=0; i<T->nrows; i++) {
						word *t0 = T->x[0]->rows[i];
						word *t1 = T->x[1]->rows[i];
						word *t2 = T->x[2]->rows[i];
						word *t3 = T->x[3]->rows[i];
						word *t4 = T->x[4]->rows[i];
						word *t5 = T->x[5]->rows[i];
						word *t6 = T->x[6]->rows[i];
						word *t7 = T->x[7]->rows[i];
						const word const *f  = F->x->rows[i];

						/* bulk of work */
						for(j=0, j2=0; j+8 < F->x->width; j+=8,j2++) {
							t0[j2] |= word_slice_64_08(f[j+0]<<7 & x80808080)>>56 | word_slice_64_08(f[j+1]<<7 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<7 & x80808080)>>40 | word_slice_64_08(f[j+3]<<7 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<7 & x80808080)>>24 | word_slice_64_08(f[j+5]<<7 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<7 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<7 & x80808080)>> 0;

							t1[j2] |= word_slice_64_08(f[j+0]<<6 & x80808080)>>56 | word_slice_64_08(f[j+1]<<6 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<6 & x80808080)>>40 | word_slice_64_08(f[j+3]<<6 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<6 & x80808080)>>24 | word_slice_64_08(f[j+5]<<6 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<6 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<6 & x80808080)>> 0;

							t2[j2] |= word_slice_64_08(f[j+0]<<5 & x80808080)>>56 | word_slice_64_08(f[j+1]<<5 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<5 & x80808080)>>40 | word_slice_64_08(f[j+3]<<5 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<5 & x80808080)>>24 | word_slice_64_08(f[j+5]<<5 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<5 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<5 & x80808080)>> 0;

							t3[j2] |= word_slice_64_08(f[j+0]<<4 & x80808080)>>56 | word_slice_64_08(f[j+1]<<4 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<4 & x80808080)>>40 | word_slice_64_08(f[j+3]<<4 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<4 & x80808080)>>24 | word_slice_64_08(f[j+5]<<4 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<4 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<4 & x80808080)>> 0;

							t4[j2] |= word_slice_64_08(f[j+0]<<3 & x80808080)>>56 | word_slice_64_08(f[j+1]<<3 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<3 & x80808080)>>40 | word_slice_64_08(f[j+3]<<3 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<3 & x80808080)>>24 | word_slice_64_08(f[j+5]<<3 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<3 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<3 & x80808080)>> 0;

							t5[j2] |= word_slice_64_08(f[j+0]<<2 & x80808080)>>56 | word_slice_64_08(f[j+1]<<2 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<2 & x80808080)>>40 | word_slice_64_08(f[j+3]<<2 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<2 & x80808080)>>24 | word_slice_64_08(f[j+5]<<2 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<2 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<2 & x80808080)>> 0;

							t6[j2] |= word_slice_64_08(f[j+0]<<1 & x80808080)>>56 | word_slice_64_08(f[j+1]<<1 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<1 & x80808080)>>40 | word_slice_64_08(f[j+3]<<1 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<1 & x80808080)>>24 | word_slice_64_08(f[j+5]<<1 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<1 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<1 & x80808080)>> 0;

							t7[j2] |= word_slice_64_08(f[j+0]<<0 & x80808080)>>56 | word_slice_64_08(f[j+1]<<0 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<0 & x80808080)>>40 | word_slice_64_08(f[j+3]<<0 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<0 & x80808080)>>24 | word_slice_64_08(f[j+5]<<0 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<0 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<0 & x80808080)>> 0;
						}
						r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0;
						switch(F->x->width - j) {
							case 8:
								r0 |= word_slice_64_08(f[j+7]<<7 & x80808080)>> 0;
								r1 |= word_slice_64_08(f[j+7]<<6 & x80808080)>> 0;
								r2 |= word_slice_64_08(f[j+7]<<5 & x80808080)>> 0;
								r3 |= word_slice_64_08(f[j+7]<<4 & x80808080)>> 0;
								r4 |= word_slice_64_08(f[j+7]<<3 & x80808080)>> 0;
								r5 |= word_slice_64_08(f[j+7]<<2 & x80808080)>> 0;
								r6 |= word_slice_64_08(f[j+7]<<1 & x80808080)>> 0;
								r7 |= word_slice_64_08(f[j+7]<<0 & x80808080)>> 0;
							case 7:
								r0 |= word_slice_64_08(f[j+6]<<7 & x80808080)>> 8;
								r1 |= word_slice_64_08(f[j+6]<<6 & x80808080)>> 8;
								r2 |= word_slice_64_08(f[j+6]<<5 & x80808080)>> 8;
								r3 |= word_slice_64_08(f[j+6]<<4 & x80808080)>> 8;
								r4 |= word_slice_64_08(f[j+6]<<3 & x80808080)>> 8;
								r5 |= word_slice_64_08(f[j+6]<<2 & x80808080)>> 8;
								r6 |= word_slice_64_08(f[j+6]<<1 & x80808080)>> 8;
								r7 |= word_slice_64_08(f[j+6]<<0 & x80808080)>> 8;
							case 6:
								r0 |= word_slice_64_08(f[j+5]<<7 & x80808080)>>16;
								r1 |= word_slice_64_08(f[j+5]<<6 & x80808080)>>16;
								r2 |= word_slice_64_08(f[j+5]<<5 & x80808080)>>16;
								r3 |= word_slice_64_08(f[j+5]<<4 & x80808080)>>16;
								r4 |= word_slice_64_08(f[j+5]<<3 & x80808080)>>16;
								r5 |= word_slice_64_08(f[j+5]<<2 & x80808080)>>16;
								r6 |= word_slice_64_08(f[j+5]<<1 & x80808080)>>16;
								r7 |= word_slice_64_08(f[j+5]<<0 & x80808080)>>16;
							case 5:
								r0 |= word_slice_64_08(f[j+4]<<7 & x80808080)>>24;
								r1 |= word_slice_64_08(f[j+4]<<6 & x80808080)>>24;
								r2 |= word_slice_64_08(f[j+4]<<5 & x80808080)>>24;
								r3 |= word_slice_64_08(f[j+4]<<4 & x80808080)>>24;
								r4 |= word_slice_64_08(f[j+4]<<3 & x80808080)>>24;
								r5 |= word_slice_64_08(f[j+4]<<2 & x80808080)>>24;
								r6 |= word_slice_64_08(f[j+4]<<1 & x80808080)>>24;
								r7 |= word_slice_64_08(f[j+4]<<0 & x80808080)>>24;
							case 4:
								r0 |= word_slice_64_08(f[j+3]<<7 & x80808080)>>32;
								r1 |= word_slice_64_08(f[j+3]<<6 & x80808080)>>32;
								r2 |= word_slice_64_08(f[j+3]<<5 & x80808080)>>32;
								r3 |= word_slice_64_08(f[j+3]<<4 & x80808080)>>32;
								r4 |= word_slice_64_08(f[j+3]<<3 & x80808080)>>32;
								r5 |= word_slice_64_08(f[j+3]<<2 & x80808080)>>32;
								r6 |= word_slice_64_08(f[j+3]<<1 & x80808080)>>32;
								r7 |= word_slice_64_08(f[j+3]<<0 & x80808080)>>32;
							case 3:
								r0 |= word_slice_64_08(f[j+2]<<7 & x80808080)>>40;
								r1 |= word_slice_64_08(f[j+2]<<6 & x80808080)>>40;
								r2 |= word_slice_64_08(f[j+2]<<5 & x80808080)>>40;
								r3 |= word_slice_64_08(f[j+2]<<4 & x80808080)>>40;
								r4 |= word_slice_64_08(f[j+2]<<3 & x80808080)>>40;
								r5 |= word_slice_64_08(f[j+2]<<2 & x80808080)>>40;
								r6 |= word_slice_64_08(f[j+2]<<1 & x80808080)>>40;
								r7 |= word_slice_64_08(f[j+2]<<0 & x80808080)>>40;
							case 2:
								r0 |= word_slice_64_08(f[j+1]<<7 & x80808080)>>48;
								r1 |= word_slice_64_08(f[j+1]<<6 & x80808080)>>48;
								r2 |= word_slice_64_08(f[j+1]<<5 & x80808080)>>48;
								r3 |= word_slice_64_08(f[j+1]<<4 & x80808080)>>48;
								r4 |= word_slice_64_08(f[j+1]<<3 & x80808080)>>48;
								r5 |= word_slice_64_08(f[j+1]<<2 & x80808080)>>48;
								r6 |= word_slice_64_08(f[j+1]<<1 & x80808080)>>48;
								r7 |= word_slice_64_08(f[j+1]<<0 & x80808080)>>48;
							case 1:
								r0 |= word_slice_64_08(f[j+0]<<7 & x80808080)>>56;
								r1 |= word_slice_64_08(f[j+0]<<6 & x80808080)>>56;
								r2 |= word_slice_64_08(f[j+0]<<5 & x80808080)>>56;
								r3 |= word_slice_64_08(f[j+0]<<4 & x80808080)>>56;
								r4 |= word_slice_64_08(f[j+0]<<3 & x80808080)>>56;
								r5 |= word_slice_64_08(f[j+0]<<2 & x80808080)>>56;
								r6 |= word_slice_64_08(f[j+0]<<1 & x80808080)>>56;
								r7 |= word_slice_64_08(f[j+0]<<0 & x80808080)>>56;
								break;
							default:
								m4ri_die("impossible");
						}
						t0[j2] |= r0 & bitmask_end;
						t1[j2] |= r1 & bitmask_end;
						t2[j2] |= r2 & bitmask_end;
						t3[j2] |= r3 & bitmask_end;
						t4[j2] |= r4 & bitmask_end;
						t5[j2] |= r5 & bitmask_end;
						t6[j2] |= r6 & bitmask_end;
						t7[j2] |= r7 & bitmask_end;
					}
				}
				break;

		case 7: {
					for(size_t i=0; i<T->nrows; i++) {
						word *t0 = T->x[0]->rows[i];
						word *t1 = T->x[1]->rows[i];
						word *t2 = T->x[2]->rows[i];
						word *t3 = T->x[3]->rows[i];
						word *t4 = T->x[4]->rows[i];
						word *t5 = T->x[5]->rows[i];
						word *t6 = T->x[6]->rows[i];
						const word const *f  = F->x->rows[i];

						/* bulk of work */
						for(j=0, j2=0; j+8 < F->x->width; j+=8,j2++) {
							t0[j2] |= word_slice_64_08(f[j+0]<<7 & x80808080)>>56 | word_slice_64_08(f[j+1]<<7 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<7 & x80808080)>>40 | word_slice_64_08(f[j+3]<<7 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<7 & x80808080)>>24 | word_slice_64_08(f[j+5]<<7 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<7 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<7 & x80808080)>> 0;

							t1[j2] |= word_slice_64_08(f[j+0]<<6 & x80808080)>>56 | word_slice_64_08(f[j+1]<<6 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<6 & x80808080)>>40 | word_slice_64_08(f[j+3]<<6 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<6 & x80808080)>>24 | word_slice_64_08(f[j+5]<<6 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<6 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<6 & x80808080)>> 0;

							t2[j2] |= word_slice_64_08(f[j+0]<<5 & x80808080)>>56 | word_slice_64_08(f[j+1]<<5 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<5 & x80808080)>>40 | word_slice_64_08(f[j+3]<<5 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<5 & x80808080)>>24 | word_slice_64_08(f[j+5]<<5 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<5 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<5 & x80808080)>> 0;

							t3[j2] |= word_slice_64_08(f[j+0]<<4 & x80808080)>>56 | word_slice_64_08(f[j+1]<<4 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<4 & x80808080)>>40 | word_slice_64_08(f[j+3]<<4 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<4 & x80808080)>>24 | word_slice_64_08(f[j+5]<<4 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<4 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<4 & x80808080)>> 0;

							t4[j2] |= word_slice_64_08(f[j+0]<<3 & x80808080)>>56 | word_slice_64_08(f[j+1]<<3 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<3 & x80808080)>>40 | word_slice_64_08(f[j+3]<<3 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<3 & x80808080)>>24 | word_slice_64_08(f[j+5]<<3 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<3 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<3 & x80808080)>> 0;

							t5[j2] |= word_slice_64_08(f[j+0]<<2 & x80808080)>>56 | word_slice_64_08(f[j+1]<<2 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<2 & x80808080)>>40 | word_slice_64_08(f[j+3]<<2 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<2 & x80808080)>>24 | word_slice_64_08(f[j+5]<<2 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<2 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<2 & x80808080)>> 0;

							t6[j2] |= word_slice_64_08(f[j+0]<<1 & x80808080)>>56 | word_slice_64_08(f[j+1]<<1 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<1 & x80808080)>>40 | word_slice_64_08(f[j+3]<<1 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<1 & x80808080)>>24 | word_slice_64_08(f[j+5]<<1 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<1 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<1 & x80808080)>> 0;
						}
						r0 = r1 = r2 = r3 = r4 = r5 = r6 = 0;
						switch(F->x->width - j) {
							case 8:
								r0 |= word_slice_64_08(f[j+7]<<7 & x80808080)>> 0;
								r1 |= word_slice_64_08(f[j+7]<<6 & x80808080)>> 0;
								r2 |= word_slice_64_08(f[j+7]<<5 & x80808080)>> 0;
								r3 |= word_slice_64_08(f[j+7]<<4 & x80808080)>> 0;
								r4 |= word_slice_64_08(f[j+7]<<3 & x80808080)>> 0;
								r5 |= word_slice_64_08(f[j+7]<<2 & x80808080)>> 0;
								r6 |= word_slice_64_08(f[j+7]<<1 & x80808080)>> 0;
							case 7:
								r0 |= word_slice_64_08(f[j+6]<<7 & x80808080)>> 8;
								r1 |= word_slice_64_08(f[j+6]<<6 & x80808080)>> 8;
								r2 |= word_slice_64_08(f[j+6]<<5 & x80808080)>> 8;
								r3 |= word_slice_64_08(f[j+6]<<4 & x80808080)>> 8;
								r4 |= word_slice_64_08(f[j+6]<<3 & x80808080)>> 8;
								r5 |= word_slice_64_08(f[j+6]<<2 & x80808080)>> 8;
								r6 |= word_slice_64_08(f[j+6]<<1 & x80808080)>> 8;
							case 6:
								r0 |= word_slice_64_08(f[j+5]<<7 & x80808080)>>16;
								r1 |= word_slice_64_08(f[j+5]<<6 & x80808080)>>16;
								r2 |= word_slice_64_08(f[j+5]<<5 & x80808080)>>16;
								r3 |= word_slice_64_08(f[j+5]<<4 & x80808080)>>16;
								r4 |= word_slice_64_08(f[j+5]<<3 & x80808080)>>16;
								r5 |= word_slice_64_08(f[j+5]<<2 & x80808080)>>16;
								r6 |= word_slice_64_08(f[j+5]<<1 & x80808080)>>16;
							case 5:
								r0 |= word_slice_64_08(f[j+4]<<7 & x80808080)>>24;
								r1 |= word_slice_64_08(f[j+4]<<6 & x80808080)>>24;
								r2 |= word_slice_64_08(f[j+4]<<5 & x80808080)>>24;
								r3 |= word_slice_64_08(f[j+4]<<4 & x80808080)>>24;
								r4 |= word_slice_64_08(f[j+4]<<3 & x80808080)>>24;
								r5 |= word_slice_64_08(f[j+4]<<2 & x80808080)>>24;
								r6 |= word_slice_64_08(f[j+4]<<1 & x80808080)>>24;
							case 4:
								r0 |= word_slice_64_08(f[j+3]<<7 & x80808080)>>32;
								r1 |= word_slice_64_08(f[j+3]<<6 & x80808080)>>32;
								r2 |= word_slice_64_08(f[j+3]<<5 & x80808080)>>32;
								r3 |= word_slice_64_08(f[j+3]<<4 & x80808080)>>32;
								r4 |= word_slice_64_08(f[j+3]<<3 & x80808080)>>32;
								r5 |= word_slice_64_08(f[j+3]<<2 & x80808080)>>32;
								r6 |= word_slice_64_08(f[j+3]<<1 & x80808080)>>32;
							case 3:
								r0 |= word_slice_64_08(f[j+2]<<7 & x80808080)>>40;
								r1 |= word_slice_64_08(f[j+2]<<6 & x80808080)>>40;
								r2 |= word_slice_64_08(f[j+2]<<5 & x80808080)>>40;
								r3 |= word_slice_64_08(f[j+2]<<4 & x80808080)>>40;
								r4 |= word_slice_64_08(f[j+2]<<3 & x80808080)>>40;
								r5 |= word_slice_64_08(f[j+2]<<2 & x80808080)>>40;
								r6 |= word_slice_64_08(f[j+2]<<1 & x80808080)>>40;
							case 2:
								r0 |= word_slice_64_08(f[j+1]<<7 & x80808080)>>48;
								r1 |= word_slice_64_08(f[j+1]<<6 & x80808080)>>48;
								r2 |= word_slice_64_08(f[j+1]<<5 & x80808080)>>48;
								r3 |= word_slice_64_08(f[j+1]<<4 & x80808080)>>48;
								r4 |= word_slice_64_08(f[j+1]<<3 & x80808080)>>48;
								r5 |= word_slice_64_08(f[j+1]<<2 & x80808080)>>48;
								r6 |= word_slice_64_08(f[j+1]<<1 & x80808080)>>48;
							case 1:
								r0 |= word_slice_64_08(f[j+0]<<7 & x80808080)>>56;
								r1 |= word_slice_64_08(f[j+0]<<6 & x80808080)>>56;
								r2 |= word_slice_64_08(f[j+0]<<5 & x80808080)>>56;
								r3 |= word_slice_64_08(f[j+0]<<4 & x80808080)>>56;
								r4 |= word_slice_64_08(f[j+0]<<3 & x80808080)>>56;
								r5 |= word_slice_64_08(f[j+0]<<2 & x80808080)>>56;
								r6 |= word_slice_64_08(f[j+0]<<1 & x80808080)>>56;
								break;
							default:
								m4ri_die("impossible");
						}
						t0[j2] |= r0 & bitmask_end;
						t1[j2] |= r1 & bitmask_end;
						t2[j2] |= r2 & bitmask_end;
						t3[j2] |= r3 & bitmask_end;
						t4[j2] |= r4 & bitmask_end;
						t5[j2] |= r5 & bitmask_end;
						t6[j2] |= r6 & bitmask_end;
					}
				}
				break;

		case 6: {
					for(size_t i=0; i<T->nrows; i++) {
						word *t0 = T->x[0]->rows[i];
						word *t1 = T->x[1]->rows[i];
						word *t2 = T->x[2]->rows[i];
						word *t3 = T->x[3]->rows[i];
						word *t4 = T->x[4]->rows[i];
						word *t5 = T->x[5]->rows[i];
						const word const *f  = F->x->rows[i];

						/* bulk of work */
						for(j=0, j2=0; j+8 < F->x->width; j+=8,j2++) {
							t0[j2] |= word_slice_64_08(f[j+0]<<7 & x80808080)>>56 | word_slice_64_08(f[j+1]<<7 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<7 & x80808080)>>40 | word_slice_64_08(f[j+3]<<7 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<7 & x80808080)>>24 | word_slice_64_08(f[j+5]<<7 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<7 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<7 & x80808080)>> 0;

							t1[j2] |= word_slice_64_08(f[j+0]<<6 & x80808080)>>56 | word_slice_64_08(f[j+1]<<6 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<6 & x80808080)>>40 | word_slice_64_08(f[j+3]<<6 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<6 & x80808080)>>24 | word_slice_64_08(f[j+5]<<6 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<6 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<6 & x80808080)>> 0;

							t2[j2] |= word_slice_64_08(f[j+0]<<5 & x80808080)>>56 | word_slice_64_08(f[j+1]<<5 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<5 & x80808080)>>40 | word_slice_64_08(f[j+3]<<5 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<5 & x80808080)>>24 | word_slice_64_08(f[j+5]<<5 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<5 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<5 & x80808080)>> 0;

							t3[j2] |= word_slice_64_08(f[j+0]<<4 & x80808080)>>56 | word_slice_64_08(f[j+1]<<4 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<4 & x80808080)>>40 | word_slice_64_08(f[j+3]<<4 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<4 & x80808080)>>24 | word_slice_64_08(f[j+5]<<4 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<4 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<4 & x80808080)>> 0;

							t4[j2] |= word_slice_64_08(f[j+0]<<3 & x80808080)>>56 | word_slice_64_08(f[j+1]<<3 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<3 & x80808080)>>40 | word_slice_64_08(f[j+3]<<3 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<3 & x80808080)>>24 | word_slice_64_08(f[j+5]<<3 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<3 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<3 & x80808080)>> 0;

							t5[j2] |= word_slice_64_08(f[j+0]<<2 & x80808080)>>56 | word_slice_64_08(f[j+1]<<2 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<2 & x80808080)>>40 | word_slice_64_08(f[j+3]<<2 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<2 & x80808080)>>24 | word_slice_64_08(f[j+5]<<2 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<2 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<2 & x80808080)>> 0;
						}
						r0 = r1 = r2 = r3 = r4 = r5 = 0;
						switch(F->x->width - j) {
							case 8:
								r0 |= word_slice_64_08(f[j+7]<<7 & x80808080)>> 0;
								r1 |= word_slice_64_08(f[j+7]<<6 & x80808080)>> 0;
								r2 |= word_slice_64_08(f[j+7]<<5 & x80808080)>> 0;
								r3 |= word_slice_64_08(f[j+7]<<4 & x80808080)>> 0;
								r4 |= word_slice_64_08(f[j+7]<<3 & x80808080)>> 0;
								r5 |= word_slice_64_08(f[j+7]<<2 & x80808080)>> 0;
							case 7:
								r0 |= word_slice_64_08(f[j+6]<<7 & x80808080)>> 8;
								r1 |= word_slice_64_08(f[j+6]<<6 & x80808080)>> 8;
								r2 |= word_slice_64_08(f[j+6]<<5 & x80808080)>> 8;
								r3 |= word_slice_64_08(f[j+6]<<4 & x80808080)>> 8;
								r4 |= word_slice_64_08(f[j+6]<<3 & x80808080)>> 8;
								r5 |= word_slice_64_08(f[j+6]<<2 & x80808080)>> 8;
							case 6:
								r0 |= word_slice_64_08(f[j+5]<<7 & x80808080)>>16;
								r1 |= word_slice_64_08(f[j+5]<<6 & x80808080)>>16;
								r2 |= word_slice_64_08(f[j+5]<<5 & x80808080)>>16;
								r3 |= word_slice_64_08(f[j+5]<<4 & x80808080)>>16;
								r4 |= word_slice_64_08(f[j+5]<<3 & x80808080)>>16;
								r5 |= word_slice_64_08(f[j+5]<<2 & x80808080)>>16;
							case 5:
								r0 |= word_slice_64_08(f[j+4]<<7 & x80808080)>>24;
								r1 |= word_slice_64_08(f[j+4]<<6 & x80808080)>>24;
								r2 |= word_slice_64_08(f[j+4]<<5 & x80808080)>>24;
								r3 |= word_slice_64_08(f[j+4]<<4 & x80808080)>>24;
								r4 |= word_slice_64_08(f[j+4]<<3 & x80808080)>>24;
								r5 |= word_slice_64_08(f[j+4]<<2 & x80808080)>>24;
							case 4:
								r0 |= word_slice_64_08(f[j+3]<<7 & x80808080)>>32;
								r1 |= word_slice_64_08(f[j+3]<<6 & x80808080)>>32;
								r2 |= word_slice_64_08(f[j+3]<<5 & x80808080)>>32;
								r3 |= word_slice_64_08(f[j+3]<<4 & x80808080)>>32;
								r4 |= word_slice_64_08(f[j+3]<<3 & x80808080)>>32;
								r5 |= word_slice_64_08(f[j+3]<<2 & x80808080)>>32;
							case 3:
								r0 |= word_slice_64_08(f[j+2]<<7 & x80808080)>>40;
								r1 |= word_slice_64_08(f[j+2]<<6 & x80808080)>>40;
								r2 |= word_slice_64_08(f[j+2]<<5 & x80808080)>>40;
								r3 |= word_slice_64_08(f[j+2]<<4 & x80808080)>>40;
								r4 |= word_slice_64_08(f[j+2]<<3 & x80808080)>>40;
								r5 |= word_slice_64_08(f[j+2]<<2 & x80808080)>>40;
							case 2:
								r0 |= word_slice_64_08(f[j+1]<<7 & x80808080)>>48;
								r1 |= word_slice_64_08(f[j+1]<<6 & x80808080)>>48;
								r2 |= word_slice_64_08(f[j+1]<<5 & x80808080)>>48;
								r3 |= word_slice_64_08(f[j+1]<<4 & x80808080)>>48;
								r4 |= word_slice_64_08(f[j+1]<<3 & x80808080)>>48;
								r5 |= word_slice_64_08(f[j+1]<<2 & x80808080)>>48;
							case 1:
								r0 |= word_slice_64_08(f[j+0]<<7 & x80808080)>>56;
								r1 |= word_slice_64_08(f[j+0]<<6 & x80808080)>>56;
								r2 |= word_slice_64_08(f[j+0]<<5 & x80808080)>>56;
								r3 |= word_slice_64_08(f[j+0]<<4 & x80808080)>>56;
								r4 |= word_slice_64_08(f[j+0]<<3 & x80808080)>>56;
								r5 |= word_slice_64_08(f[j+0]<<2 & x80808080)>>56;
								break;
							default:
								m4ri_die("impossible");
						}
						t0[j2] |= r0 & bitmask_end;
						t1[j2] |= r1 & bitmask_end;
						t2[j2] |= r2 & bitmask_end;
						t3[j2] |= r3 & bitmask_end;
						t4[j2] |= r4 & bitmask_end;
						t5[j2] |= r5 & bitmask_end;
					}
				}
				break;

		case 5: {
					for(size_t i=0; i<T->nrows; i++) {
						word *t0 = T->x[0]->rows[i];
						word *t1 = T->x[1]->rows[i];
						word *t2 = T->x[2]->rows[i];
						word *t3 = T->x[3]->rows[i];
						word *t4 = T->x[4]->rows[i];
						const word const *f  = F->x->rows[i];

						/* bulk of work */
						for(j=0, j2=0; j+8 < F->x->width; j+=8,j2++) {
							t0[j2] |= word_slice_64_08(f[j+0]<<7 & x80808080)>>56 | word_slice_64_08(f[j+1]<<7 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<7 & x80808080)>>40 | word_slice_64_08(f[j+3]<<7 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<7 & x80808080)>>24 | word_slice_64_08(f[j+5]<<7 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<7 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<7 & x80808080)>> 0;

							t1[j2] |= word_slice_64_08(f[j+0]<<6 & x80808080)>>56 | word_slice_64_08(f[j+1]<<6 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<6 & x80808080)>>40 | word_slice_64_08(f[j+3]<<6 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<6 & x80808080)>>24 | word_slice_64_08(f[j+5]<<6 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<6 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<6 & x80808080)>> 0;

							t2[j2] |= word_slice_64_08(f[j+0]<<5 & x80808080)>>56 | word_slice_64_08(f[j+1]<<5 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<5 & x80808080)>>40 | word_slice_64_08(f[j+3]<<5 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<5 & x80808080)>>24 | word_slice_64_08(f[j+5]<<5 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<5 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<5 & x80808080)>> 0;

							t3[j2] |= word_slice_64_08(f[j+0]<<4 & x80808080)>>56 | word_slice_64_08(f[j+1]<<4 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<4 & x80808080)>>40 | word_slice_64_08(f[j+3]<<4 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<4 & x80808080)>>24 | word_slice_64_08(f[j+5]<<4 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<4 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<4 & x80808080)>> 0;

							t4[j2] |= word_slice_64_08(f[j+0]<<3 & x80808080)>>56 | word_slice_64_08(f[j+1]<<3 & x80808080)>>48 \
									  |       word_slice_64_08(f[j+2]<<3 & x80808080)>>40 | word_slice_64_08(f[j+3]<<3 & x80808080)>>32 \
									  |       word_slice_64_08(f[j+4]<<3 & x80808080)>>24 | word_slice_64_08(f[j+5]<<3 & x80808080)>>16 \
									  |       word_slice_64_08(f[j+6]<<3 & x80808080)>> 8 | word_slice_64_08(f[j+7]<<3 & x80808080)>> 0;
						}
						r0 = r1 = r2 = r3 = r4 = 0;
						switch(F->x->width - j) {
							case 8:
								r0 |= word_slice_64_08(f[j+7]<<7 & x80808080)>> 0;
								r1 |= word_slice_64_08(f[j+7]<<6 & x80808080)>> 0;
								r2 |= word_slice_64_08(f[j+7]<<5 & x80808080)>> 0;
								r3 |= word_slice_64_08(f[j+7]<<4 & x80808080)>> 0;
								r4 |= word_slice_64_08(f[j+7]<<3 & x80808080)>> 0;
							case 7:
								r0 |= word_slice_64_08(f[j+6]<<7 & x80808080)>> 8;
								r1 |= word_slice_64_08(f[j+6]<<6 & x80808080)>> 8;
								r2 |= word_slice_64_08(f[j+6]<<5 & x80808080)>> 8;
								r3 |= word_slice_64_08(f[j+6]<<4 & x80808080)>> 8;
								r4 |= word_slice_64_08(f[j+6]<<3 & x80808080)>> 8;
							case 6:
								r0 |= word_slice_64_08(f[j+5]<<7 & x80808080)>>16;
								r1 |= word_slice_64_08(f[j+5]<<6 & x80808080)>>16;
								r2 |= word_slice_64_08(f[j+5]<<5 & x80808080)>>16;
								r3 |= word_slice_64_08(f[j+5]<<4 & x80808080)>>16;
								r4 |= word_slice_64_08(f[j+5]<<3 & x80808080)>>16;
							case 5:
								r0 |= word_slice_64_08(f[j+4]<<7 & x80808080)>>24;
								r1 |= word_slice_64_08(f[j+4]<<6 & x80808080)>>24;
								r2 |= word_slice_64_08(f[j+4]<<5 & x80808080)>>24;
								r3 |= word_slice_64_08(f[j+4]<<4 & x80808080)>>24;
								r4 |= word_slice_64_08(f[j+4]<<3 & x80808080)>>24;
							case 4:
								r0 |= word_slice_64_08(f[j+3]<<7 & x80808080)>>32;
								r1 |= word_slice_64_08(f[j+3]<<6 & x80808080)>>32;
								r2 |= word_slice_64_08(f[j+3]<<5 & x80808080)>>32;
								r3 |= word_slice_64_08(f[j+3]<<4 & x80808080)>>32;
								r4 |= word_slice_64_08(f[j+3]<<3 & x80808080)>>32;
							case 3:
								r0 |= word_slice_64_08(f[j+2]<<7 & x80808080)>>40;
								r1 |= word_slice_64_08(f[j+2]<<6 & x80808080)>>40;
								r2 |= word_slice_64_08(f[j+2]<<5 & x80808080)>>40;
								r3 |= word_slice_64_08(f[j+2]<<4 & x80808080)>>40;
								r4 |= word_slice_64_08(f[j+2]<<3 & x80808080)>>40;
							case 2:
								r0 |= word_slice_64_08(f[j+1]<<7 & x80808080)>>48;
								r1 |= word_slice_64_08(f[j+1]<<6 & x80808080)>>48;
								r2 |= word_slice_64_08(f[j+1]<<5 & x80808080)>>48;
								r3 |= word_slice_64_08(f[j+1]<<4 & x80808080)>>48;
								r4 |= word_slice_64_08(f[j+1]<<3 & x80808080)>>48;
							case 1:
								r0 |= word_slice_64_08(f[j+0]<<7 & x80808080)>>56;
								r1 |= word_slice_64_08(f[j+0]<<6 & x80808080)>>56;
								r2 |= word_slice_64_08(f[j+0]<<5 & x80808080)>>56;
								r3 |= word_slice_64_08(f[j+0]<<4 & x80808080)>>56;
								r4 |= word_slice_64_08(f[j+0]<<3 & x80808080)>>56;
								break;
							default:
								m4ri_die("impossible");
						}
						t0[j2] |= r0 & bitmask_end;
						t1[j2] |= r1 & bitmask_end;
						t2[j2] |= r2 & bitmask_end;
						t3[j2] |= r3 & bitmask_end;
						t4[j2] |= r4 & bitmask_end;
					}
				}
				break;

		default:
				m4ri_die("impossible\n");
	}
	return T;
}
