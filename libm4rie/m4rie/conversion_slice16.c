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

static inline word word_slice_64_16(word a) {
	a = (a & xffff0000) | (a & xffff0000>>16)<<15;
	a = (a & xffffffff) | (a & xffffffff>>32)<<30;
	return a;
}

/* we define these things to keep code compact below. */

#define word_slice_64_16_combine_bulk(T, Ti, F, Fi, shift)  \
	T[Ti] |= word_slice_64_16(F[Fi+ 0]<<shift & x80008000)>>60 | word_slice_64_16(F[Fi+ 1]<<shift & x80008000)>>56 \
|       word_slice_64_16(F[Fi+ 2]<<shift & x80008000)>>52 | word_slice_64_16(F[Fi+ 3]<<shift & x80008000)>>48 \
|       word_slice_64_16(F[Fi+ 4]<<shift & x80008000)>>44 | word_slice_64_16(F[Fi+ 5]<<shift & x80008000)>>40 \
|       word_slice_64_16(F[Fi+ 6]<<shift & x80008000)>>36 | word_slice_64_16(F[Fi+ 7]<<shift & x80008000)>>32 \
|       word_slice_64_16(F[Fi+ 8]<<shift & x80008000)>>28 | word_slice_64_16(F[Fi+ 9]<<shift & x80008000)>>24 \
|       word_slice_64_16(F[Fi+10]<<shift & x80008000)>>20 | word_slice_64_16(F[Fi+11]<<shift & x80008000)>>16 \
|       word_slice_64_16(F[Fi+12]<<shift & x80008000)>>12 | word_slice_64_16(F[Fi+13]<<shift & x80008000)>> 8 \
|       word_slice_64_16(F[Fi+14]<<shift & x80008000)>> 4 | word_slice_64_16(F[Fi+15]<<shift & x80008000)>> 0;

#define word_slice_64_16_slice_rest(F, Fi, shift)         \
	r0 |= word_slice_64_16(F[Fi]<<15 & x80008000)>> shift;         \
r1 |= word_slice_64_16(F[Fi]<<14 & x80008000)>> shift;         \
r2 |= word_slice_64_16(F[Fi]<<13 & x80008000)>> shift;         \
r3 |= word_slice_64_16(F[Fi]<<12 & x80008000)>> shift;         \
r4 |= word_slice_64_16(F[Fi]<<11 & x80008000)>> shift;         \
r5 |= word_slice_64_16(F[Fi]<<10 & x80008000)>> shift;         \
r6 |= word_slice_64_16(F[Fi]<< 9 & x80008000)>> shift;         \
r7 |= word_slice_64_16(F[Fi]<< 8 & x80008000)>> shift;

mzd_slice_t *_mzed_slice16(mzd_slice_t *T, const mzed_t *F) {
	assert(T && (8 < T->depth && T->depth <= 16));
	size_t j, j2 = 0;
	register word r0,r1,r2,r3,r4,r5,r6,r7 = 0;

	const word bitmask_end = T->x[0]->high_bitmask;

	if (mzed_is_zero(F))
		return T;

	/* we do multiple runs over T to make the code more compact, we start by doing the first eight
	   bits */

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
		for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
			word_slice_64_16_combine_bulk(t0, j2, f, j, 15);
			word_slice_64_16_combine_bulk(t1, j2, f, j, 14);
			word_slice_64_16_combine_bulk(t2, j2, f, j, 13);
			word_slice_64_16_combine_bulk(t3, j2, f, j, 12);
			word_slice_64_16_combine_bulk(t4, j2, f, j, 11);
			word_slice_64_16_combine_bulk(t5, j2, f, j, 10);
			word_slice_64_16_combine_bulk(t6, j2, f, j,  9);
			word_slice_64_16_combine_bulk(t7, j2, f, j,  8);
		}
		r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0;
		switch(F->x->width - j) {
			case 16: word_slice_64_16_slice_rest(f, j+15,  0);
			case 15: word_slice_64_16_slice_rest(f, j+14,  4);
			case 14: word_slice_64_16_slice_rest(f, j+13,  8);
			case 13: word_slice_64_16_slice_rest(f, j+12, 12);
			case 12: word_slice_64_16_slice_rest(f, j+11, 16);
			case 11: word_slice_64_16_slice_rest(f, j+10, 20);
			case 10: word_slice_64_16_slice_rest(f, j+ 9, 24);
			case  9: word_slice_64_16_slice_rest(f, j+ 8, 28);
			case  8: word_slice_64_16_slice_rest(f, j+ 7, 32);
			case  7: word_slice_64_16_slice_rest(f, j+ 6, 36);
			case  6: word_slice_64_16_slice_rest(f, j+ 5, 40);
			case  5: word_slice_64_16_slice_rest(f, j+ 4, 44);
			case  4: word_slice_64_16_slice_rest(f, j+ 3, 48);
			case  3: word_slice_64_16_slice_rest(f, j+ 2, 52);
			case  2: word_slice_64_16_slice_rest(f, j+ 1, 56);
			case  1: word_slice_64_16_slice_rest(f, j+ 0, 60);
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
	if(T->depth >= 12) {
		for(size_t i=0; i<T->nrows; i++) {
			word *t0 = T->x[ 8]->rows[i];
			word *t1 = T->x[ 9]->rows[i];
			word *t2 = T->x[10]->rows[i];
			word *t3 = T->x[11]->rows[i];
			const word const *f  = F->x->rows[i];

			/* bulk of work */
			for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
				word_slice_64_16_combine_bulk(t0, j2, f, j,  7);
				word_slice_64_16_combine_bulk(t1, j2, f, j,  6);
				word_slice_64_16_combine_bulk(t2, j2, f, j,  5);
				word_slice_64_16_combine_bulk(t3, j2, f, j,  4);
			}
			r0 = r1 = r2 = r3 = 0;
			switch(F->x->width - j) {
				case 16: r0 |= word_slice_64_16(f[j+15]<< 7 & x80008000)>>  0; r1 |= word_slice_64_16(f[j+15]<< 6 & x80008000)>>  0; r2 |= word_slice_64_16(f[j+15]<< 5 & x80008000)>>  0; r3 |= word_slice_64_16(f[j+15]<< 4 & x80008000)>>  0;
				case 15: r0 |= word_slice_64_16(f[j+14]<< 7 & x80008000)>>  4; r1 |= word_slice_64_16(f[j+14]<< 6 & x80008000)>>  4; r2 |= word_slice_64_16(f[j+14]<< 5 & x80008000)>>  4; r3 |= word_slice_64_16(f[j+14]<< 4 & x80008000)>>  4;
				case 14: r0 |= word_slice_64_16(f[j+13]<< 7 & x80008000)>>  8; r1 |= word_slice_64_16(f[j+13]<< 6 & x80008000)>>  8; r2 |= word_slice_64_16(f[j+13]<< 5 & x80008000)>>  8; r3 |= word_slice_64_16(f[j+13]<< 4 & x80008000)>>  8;
				case 13: r0 |= word_slice_64_16(f[j+12]<< 7 & x80008000)>> 12; r1 |= word_slice_64_16(f[j+12]<< 6 & x80008000)>> 12; r2 |= word_slice_64_16(f[j+12]<< 5 & x80008000)>> 12; r3 |= word_slice_64_16(f[j+12]<< 4 & x80008000)>> 12;
				case 12: r0 |= word_slice_64_16(f[j+11]<< 7 & x80008000)>> 16; r1 |= word_slice_64_16(f[j+11]<< 6 & x80008000)>> 16; r2 |= word_slice_64_16(f[j+11]<< 5 & x80008000)>> 16; r3 |= word_slice_64_16(f[j+11]<< 4 & x80008000)>> 16;
				case 11: r0 |= word_slice_64_16(f[j+10]<< 7 & x80008000)>> 20; r1 |= word_slice_64_16(f[j+10]<< 6 & x80008000)>> 20; r2 |= word_slice_64_16(f[j+10]<< 5 & x80008000)>> 20; r3 |= word_slice_64_16(f[j+10]<< 4 & x80008000)>> 20;
				case 10: r0 |= word_slice_64_16(f[j+ 9]<< 7 & x80008000)>> 24; r1 |= word_slice_64_16(f[j+ 9]<< 6 & x80008000)>> 24; r2 |= word_slice_64_16(f[j+ 9]<< 5 & x80008000)>> 24; r3 |= word_slice_64_16(f[j+ 9]<< 4 & x80008000)>> 24;
				case  9: r0 |= word_slice_64_16(f[j+ 8]<< 7 & x80008000)>> 28; r1 |= word_slice_64_16(f[j+ 8]<< 6 & x80008000)>> 28; r2 |= word_slice_64_16(f[j+ 8]<< 5 & x80008000)>> 28; r3 |= word_slice_64_16(f[j+ 8]<< 4 & x80008000)>> 28;
				case  8: r0 |= word_slice_64_16(f[j+ 7]<< 7 & x80008000)>> 32; r1 |= word_slice_64_16(f[j+ 7]<< 6 & x80008000)>> 32; r2 |= word_slice_64_16(f[j+ 7]<< 5 & x80008000)>> 32; r3 |= word_slice_64_16(f[j+ 7]<< 4 & x80008000)>> 32;
				case  7: r0 |= word_slice_64_16(f[j+ 6]<< 7 & x80008000)>> 36; r1 |= word_slice_64_16(f[j+ 6]<< 6 & x80008000)>> 36; r2 |= word_slice_64_16(f[j+ 6]<< 5 & x80008000)>> 36; r3 |= word_slice_64_16(f[j+ 6]<< 4 & x80008000)>> 36;
				case  6: r0 |= word_slice_64_16(f[j+ 5]<< 7 & x80008000)>> 40; r1 |= word_slice_64_16(f[j+ 5]<< 6 & x80008000)>> 40; r2 |= word_slice_64_16(f[j+ 5]<< 5 & x80008000)>> 40; r3 |= word_slice_64_16(f[j+ 5]<< 4 & x80008000)>> 40;
				case  5: r0 |= word_slice_64_16(f[j+ 4]<< 7 & x80008000)>> 44; r1 |= word_slice_64_16(f[j+ 4]<< 6 & x80008000)>> 44; r2 |= word_slice_64_16(f[j+ 4]<< 5 & x80008000)>> 44; r3 |= word_slice_64_16(f[j+ 4]<< 4 & x80008000)>> 44;
				case  4: r0 |= word_slice_64_16(f[j+ 3]<< 7 & x80008000)>> 48; r1 |= word_slice_64_16(f[j+ 3]<< 6 & x80008000)>> 48; r2 |= word_slice_64_16(f[j+ 3]<< 5 & x80008000)>> 48; r3 |= word_slice_64_16(f[j+ 3]<< 4 & x80008000)>> 48;
				case  3: r0 |= word_slice_64_16(f[j+ 2]<< 7 & x80008000)>> 52; r1 |= word_slice_64_16(f[j+ 2]<< 6 & x80008000)>> 52; r2 |= word_slice_64_16(f[j+ 2]<< 5 & x80008000)>> 52; r3 |= word_slice_64_16(f[j+ 2]<< 4 & x80008000)>> 52;
				case  2: r0 |= word_slice_64_16(f[j+ 1]<< 7 & x80008000)>> 56; r1 |= word_slice_64_16(f[j+ 1]<< 6 & x80008000)>> 56; r2 |= word_slice_64_16(f[j+ 1]<< 5 & x80008000)>> 56; r3 |= word_slice_64_16(f[j+ 1]<< 4 & x80008000)>> 56;
				case  1: r0 |= word_slice_64_16(f[j+ 0]<< 7 & x80008000)>> 60; r1 |= word_slice_64_16(f[j+ 0]<< 6 & x80008000)>> 60; r2 |= word_slice_64_16(f[j+ 0]<< 5 & x80008000)>> 60; r3 |= word_slice_64_16(f[j+ 0]<< 4 & x80008000)>> 60;
						 break;
				default:
						 m4ri_die("impossible");
			}
			t0[j2] |= r0 & bitmask_end;
			t1[j2] |= r1 & bitmask_end;
			t2[j2] |= r2 & bitmask_end;
			t3[j2] |= r3 & bitmask_end;
		}

		switch(T->depth) {
			case 16: {
						 for(size_t i=0; i<T->nrows; i++) {
							 word *t0 = T->x[12]->rows[i];
							 word *t1 = T->x[13]->rows[i];
							 word *t2 = T->x[14]->rows[i];
							 word *t3 = T->x[15]->rows[i];
							 const word const *f  = F->x->rows[i];

							 /* bulk of work */
							 for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
								 word_slice_64_16_combine_bulk(t0, j2, f, j,  3);
								 word_slice_64_16_combine_bulk(t1, j2, f, j,  2);
								 word_slice_64_16_combine_bulk(t2, j2, f, j,  1);
								 word_slice_64_16_combine_bulk(t3, j2, f, j,  0);
							 }
							 r0 = r1 = r2 = r3 = 0;
							 switch(F->x->width - j) {
								 case 16: r0 |= word_slice_64_16(f[j+15]<< 3 & x80008000)>>  0; r1 |= word_slice_64_16(f[j+15]<< 2 & x80008000)>>  0; r2 |= word_slice_64_16(f[j+15]<< 1 & x80008000)>>  0; r3 |= word_slice_64_16(f[j+15]<< 0 & x80008000)>>  0;
								 case 15: r0 |= word_slice_64_16(f[j+14]<< 3 & x80008000)>>  4; r1 |= word_slice_64_16(f[j+14]<< 2 & x80008000)>>  4; r2 |= word_slice_64_16(f[j+14]<< 1 & x80008000)>>  4; r3 |= word_slice_64_16(f[j+14]<< 0 & x80008000)>>  4;
								 case 14: r0 |= word_slice_64_16(f[j+13]<< 3 & x80008000)>>  8; r1 |= word_slice_64_16(f[j+13]<< 2 & x80008000)>>  8; r2 |= word_slice_64_16(f[j+13]<< 1 & x80008000)>>  8; r3 |= word_slice_64_16(f[j+13]<< 0 & x80008000)>>  8;
								 case 13: r0 |= word_slice_64_16(f[j+12]<< 3 & x80008000)>> 12; r1 |= word_slice_64_16(f[j+12]<< 2 & x80008000)>> 12; r2 |= word_slice_64_16(f[j+12]<< 1 & x80008000)>> 12; r3 |= word_slice_64_16(f[j+12]<< 0 & x80008000)>> 12;
								 case 12: r0 |= word_slice_64_16(f[j+11]<< 3 & x80008000)>> 16; r1 |= word_slice_64_16(f[j+11]<< 2 & x80008000)>> 16; r2 |= word_slice_64_16(f[j+11]<< 1 & x80008000)>> 16; r3 |= word_slice_64_16(f[j+11]<< 0 & x80008000)>> 16;
								 case 11: r0 |= word_slice_64_16(f[j+10]<< 3 & x80008000)>> 20; r1 |= word_slice_64_16(f[j+10]<< 2 & x80008000)>> 20; r2 |= word_slice_64_16(f[j+10]<< 1 & x80008000)>> 20; r3 |= word_slice_64_16(f[j+10]<< 0 & x80008000)>> 20;
								 case 10: r0 |= word_slice_64_16(f[j+ 9]<< 3 & x80008000)>> 24; r1 |= word_slice_64_16(f[j+ 9]<< 2 & x80008000)>> 24; r2 |= word_slice_64_16(f[j+ 9]<< 1 & x80008000)>> 24; r3 |= word_slice_64_16(f[j+ 9]<< 0 & x80008000)>> 24;
								 case  9: r0 |= word_slice_64_16(f[j+ 8]<< 3 & x80008000)>> 28; r1 |= word_slice_64_16(f[j+ 8]<< 2 & x80008000)>> 28; r2 |= word_slice_64_16(f[j+ 8]<< 1 & x80008000)>> 28; r3 |= word_slice_64_16(f[j+ 8]<< 0 & x80008000)>> 28;
								 case  8: r0 |= word_slice_64_16(f[j+ 7]<< 3 & x80008000)>> 32; r1 |= word_slice_64_16(f[j+ 7]<< 2 & x80008000)>> 32; r2 |= word_slice_64_16(f[j+ 7]<< 1 & x80008000)>> 32; r3 |= word_slice_64_16(f[j+ 7]<< 0 & x80008000)>> 32;
								 case  7: r0 |= word_slice_64_16(f[j+ 6]<< 3 & x80008000)>> 36; r1 |= word_slice_64_16(f[j+ 6]<< 2 & x80008000)>> 36; r2 |= word_slice_64_16(f[j+ 6]<< 1 & x80008000)>> 36; r3 |= word_slice_64_16(f[j+ 6]<< 0 & x80008000)>> 36;
								 case  6: r0 |= word_slice_64_16(f[j+ 5]<< 3 & x80008000)>> 40; r1 |= word_slice_64_16(f[j+ 5]<< 2 & x80008000)>> 40; r2 |= word_slice_64_16(f[j+ 5]<< 1 & x80008000)>> 40; r3 |= word_slice_64_16(f[j+ 5]<< 0 & x80008000)>> 40;
								 case  5: r0 |= word_slice_64_16(f[j+ 4]<< 3 & x80008000)>> 44; r1 |= word_slice_64_16(f[j+ 4]<< 2 & x80008000)>> 44; r2 |= word_slice_64_16(f[j+ 4]<< 1 & x80008000)>> 44; r3 |= word_slice_64_16(f[j+ 4]<< 0 & x80008000)>> 44;
								 case  4: r0 |= word_slice_64_16(f[j+ 3]<< 3 & x80008000)>> 48; r1 |= word_slice_64_16(f[j+ 3]<< 2 & x80008000)>> 48; r2 |= word_slice_64_16(f[j+ 3]<< 1 & x80008000)>> 48; r3 |= word_slice_64_16(f[j+ 3]<< 0 & x80008000)>> 48;
								 case  3: r0 |= word_slice_64_16(f[j+ 2]<< 3 & x80008000)>> 52; r1 |= word_slice_64_16(f[j+ 2]<< 2 & x80008000)>> 52; r2 |= word_slice_64_16(f[j+ 2]<< 1 & x80008000)>> 52; r3 |= word_slice_64_16(f[j+ 2]<< 0 & x80008000)>> 52;
								 case  2: r0 |= word_slice_64_16(f[j+ 1]<< 3 & x80008000)>> 56; r1 |= word_slice_64_16(f[j+ 1]<< 2 & x80008000)>> 56; r2 |= word_slice_64_16(f[j+ 1]<< 1 & x80008000)>> 56; r3 |= word_slice_64_16(f[j+ 1]<< 0 & x80008000)>> 56;
								 case  1: r0 |= word_slice_64_16(f[j+ 0]<< 3 & x80008000)>> 60; r1 |= word_slice_64_16(f[j+ 0]<< 2 & x80008000)>> 60; r2 |= word_slice_64_16(f[j+ 0]<< 1 & x80008000)>> 60; r3 |= word_slice_64_16(f[j+ 0]<< 0 & x80008000)>> 60;
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t0[j2] |= r0 & bitmask_end;
							 t1[j2] |= r1 & bitmask_end;
							 t2[j2] |= r2 & bitmask_end;
							 t3[j2] |= r3 & bitmask_end;
						 }
					 } break;
			case 15: {
						 for(size_t i=0; i<T->nrows; i++) {
							 word *t0 = T->x[12]->rows[i];
							 word *t1 = T->x[13]->rows[i];
							 word *t2 = T->x[14]->rows[i];
							 const word const *f  = F->x->rows[i];

							 /* bulk of work */
							 for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
								 word_slice_64_16_combine_bulk(t0, j2, f, j,  3);
								 word_slice_64_16_combine_bulk(t1, j2, f, j,  2);
								 word_slice_64_16_combine_bulk(t2, j2, f, j,  1);
							 }
							 r0 = r1 = r2 = 0;
							 switch(F->x->width - j) {
								 case 16: r0 |= word_slice_64_16(f[j+15]<< 3 & x80008000)>>  0; r1 |= word_slice_64_16(f[j+15]<< 2 & x80008000)>>  0; r2 |= word_slice_64_16(f[j+15]<< 1 & x80008000)>>  0;
								 case 15: r0 |= word_slice_64_16(f[j+14]<< 3 & x80008000)>>  4; r1 |= word_slice_64_16(f[j+14]<< 2 & x80008000)>>  4; r2 |= word_slice_64_16(f[j+14]<< 1 & x80008000)>>  4;
								 case 14: r0 |= word_slice_64_16(f[j+13]<< 3 & x80008000)>>  8; r1 |= word_slice_64_16(f[j+13]<< 2 & x80008000)>>  8; r2 |= word_slice_64_16(f[j+13]<< 1 & x80008000)>>  8;
								 case 13: r0 |= word_slice_64_16(f[j+12]<< 3 & x80008000)>> 12; r1 |= word_slice_64_16(f[j+12]<< 2 & x80008000)>> 12; r2 |= word_slice_64_16(f[j+12]<< 1 & x80008000)>> 12;
								 case 12: r0 |= word_slice_64_16(f[j+11]<< 3 & x80008000)>> 16; r1 |= word_slice_64_16(f[j+11]<< 2 & x80008000)>> 16; r2 |= word_slice_64_16(f[j+11]<< 1 & x80008000)>> 16;
								 case 11: r0 |= word_slice_64_16(f[j+10]<< 3 & x80008000)>> 20; r1 |= word_slice_64_16(f[j+10]<< 2 & x80008000)>> 20; r2 |= word_slice_64_16(f[j+10]<< 1 & x80008000)>> 20;
								 case 10: r0 |= word_slice_64_16(f[j+ 9]<< 3 & x80008000)>> 24; r1 |= word_slice_64_16(f[j+ 9]<< 2 & x80008000)>> 24; r2 |= word_slice_64_16(f[j+ 9]<< 1 & x80008000)>> 24;
								 case  9: r0 |= word_slice_64_16(f[j+ 8]<< 3 & x80008000)>> 28; r1 |= word_slice_64_16(f[j+ 8]<< 2 & x80008000)>> 28; r2 |= word_slice_64_16(f[j+ 8]<< 1 & x80008000)>> 28;
								 case  8: r0 |= word_slice_64_16(f[j+ 7]<< 3 & x80008000)>> 32; r1 |= word_slice_64_16(f[j+ 7]<< 2 & x80008000)>> 32; r2 |= word_slice_64_16(f[j+ 7]<< 1 & x80008000)>> 32;
								 case  7: r0 |= word_slice_64_16(f[j+ 6]<< 3 & x80008000)>> 36; r1 |= word_slice_64_16(f[j+ 6]<< 2 & x80008000)>> 36; r2 |= word_slice_64_16(f[j+ 6]<< 1 & x80008000)>> 36;
								 case  6: r0 |= word_slice_64_16(f[j+ 5]<< 3 & x80008000)>> 40; r1 |= word_slice_64_16(f[j+ 5]<< 2 & x80008000)>> 40; r2 |= word_slice_64_16(f[j+ 5]<< 1 & x80008000)>> 40;
								 case  5: r0 |= word_slice_64_16(f[j+ 4]<< 3 & x80008000)>> 44; r1 |= word_slice_64_16(f[j+ 4]<< 2 & x80008000)>> 44; r2 |= word_slice_64_16(f[j+ 4]<< 1 & x80008000)>> 44;
								 case  4: r0 |= word_slice_64_16(f[j+ 3]<< 3 & x80008000)>> 48; r1 |= word_slice_64_16(f[j+ 3]<< 2 & x80008000)>> 48; r2 |= word_slice_64_16(f[j+ 3]<< 1 & x80008000)>> 48;
								 case  3: r0 |= word_slice_64_16(f[j+ 2]<< 3 & x80008000)>> 52; r1 |= word_slice_64_16(f[j+ 2]<< 2 & x80008000)>> 52; r2 |= word_slice_64_16(f[j+ 2]<< 1 & x80008000)>> 52;
								 case  2: r0 |= word_slice_64_16(f[j+ 1]<< 3 & x80008000)>> 56; r1 |= word_slice_64_16(f[j+ 1]<< 2 & x80008000)>> 56; r2 |= word_slice_64_16(f[j+ 1]<< 1 & x80008000)>> 56;
								 case  1: r0 |= word_slice_64_16(f[j+ 0]<< 3 & x80008000)>> 60; r1 |= word_slice_64_16(f[j+ 0]<< 2 & x80008000)>> 60; r2 |= word_slice_64_16(f[j+ 0]<< 1 & x80008000)>> 60;
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t0[j2] |= r0 & bitmask_end;
							 t1[j2] |= r1 & bitmask_end;
							 t2[j2] |= r2 & bitmask_end;
						 }
					 } break;
			case 14: {
						 for(size_t i=0; i<T->nrows; i++) {
							 word *t0 = T->x[12]->rows[i];
							 word *t1 = T->x[13]->rows[i];
							 const word const *f  = F->x->rows[i];

							 /* bulk of work */
							 for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
								 word_slice_64_16_combine_bulk(t0, j2, f, j,  3);
								 word_slice_64_16_combine_bulk(t1, j2, f, j,  2);
							 }
							 r0 = r1 = 0;
							 switch(F->x->width - j) {
								 case 16: r0 |= word_slice_64_16(f[j+15]<< 3 & x80008000)>>  0; r1 |= word_slice_64_16(f[j+15]<< 2 & x80008000)>>  0;
								 case 15: r0 |= word_slice_64_16(f[j+14]<< 3 & x80008000)>>  4; r1 |= word_slice_64_16(f[j+14]<< 2 & x80008000)>>  4;
								 case 14: r0 |= word_slice_64_16(f[j+13]<< 3 & x80008000)>>  8; r1 |= word_slice_64_16(f[j+13]<< 2 & x80008000)>>  8;
								 case 13: r0 |= word_slice_64_16(f[j+12]<< 3 & x80008000)>> 12; r1 |= word_slice_64_16(f[j+12]<< 2 & x80008000)>> 12;
								 case 12: r0 |= word_slice_64_16(f[j+11]<< 3 & x80008000)>> 16; r1 |= word_slice_64_16(f[j+11]<< 2 & x80008000)>> 16;
								 case 11: r0 |= word_slice_64_16(f[j+10]<< 3 & x80008000)>> 20; r1 |= word_slice_64_16(f[j+10]<< 2 & x80008000)>> 20;
								 case 10: r0 |= word_slice_64_16(f[j+ 9]<< 3 & x80008000)>> 24; r1 |= word_slice_64_16(f[j+ 9]<< 2 & x80008000)>> 24;
								 case  9: r0 |= word_slice_64_16(f[j+ 8]<< 3 & x80008000)>> 28; r1 |= word_slice_64_16(f[j+ 8]<< 2 & x80008000)>> 28;
								 case  8: r0 |= word_slice_64_16(f[j+ 7]<< 3 & x80008000)>> 32; r1 |= word_slice_64_16(f[j+ 7]<< 2 & x80008000)>> 32;
								 case  7: r0 |= word_slice_64_16(f[j+ 6]<< 3 & x80008000)>> 36; r1 |= word_slice_64_16(f[j+ 6]<< 2 & x80008000)>> 36;
								 case  6: r0 |= word_slice_64_16(f[j+ 5]<< 3 & x80008000)>> 40; r1 |= word_slice_64_16(f[j+ 5]<< 2 & x80008000)>> 40;
								 case  5: r0 |= word_slice_64_16(f[j+ 4]<< 3 & x80008000)>> 44; r1 |= word_slice_64_16(f[j+ 4]<< 2 & x80008000)>> 44;
								 case  4: r0 |= word_slice_64_16(f[j+ 3]<< 3 & x80008000)>> 48; r1 |= word_slice_64_16(f[j+ 3]<< 2 & x80008000)>> 48;
								 case  3: r0 |= word_slice_64_16(f[j+ 2]<< 3 & x80008000)>> 52; r1 |= word_slice_64_16(f[j+ 2]<< 2 & x80008000)>> 52;
								 case  2: r0 |= word_slice_64_16(f[j+ 1]<< 3 & x80008000)>> 56; r1 |= word_slice_64_16(f[j+ 1]<< 2 & x80008000)>> 56;
								 case  1: r0 |= word_slice_64_16(f[j+ 0]<< 3 & x80008000)>> 60; r1 |= word_slice_64_16(f[j+ 0]<< 2 & x80008000)>> 60;
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t0[j2] |= r0 & bitmask_end;
							 t1[j2] |= r1 & bitmask_end;
						 }
					 } break;
			case 13: {
						 for(size_t i=0; i<T->nrows; i++) {
							 word *t0 = T->x[12]->rows[i];
							 const word const *f  = F->x->rows[i];

							 /* bulk of work */
							 for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
								 word_slice_64_16_combine_bulk(t0, j2, f, j,  3);
							 }
							 r0 = 0;
							 switch(F->x->width - j) {
								 case 16: r0 |= word_slice_64_16(f[j+15]<< 3 & x80008000)>>  0;
								 case 15: r0 |= word_slice_64_16(f[j+14]<< 3 & x80008000)>>  4;
								 case 14: r0 |= word_slice_64_16(f[j+13]<< 3 & x80008000)>>  8;
								 case 13: r0 |= word_slice_64_16(f[j+12]<< 3 & x80008000)>> 12;
								 case 12: r0 |= word_slice_64_16(f[j+11]<< 3 & x80008000)>> 16;
								 case 11: r0 |= word_slice_64_16(f[j+10]<< 3 & x80008000)>> 20;
								 case 10: r0 |= word_slice_64_16(f[j+ 9]<< 3 & x80008000)>> 24;
								 case  9: r0 |= word_slice_64_16(f[j+ 8]<< 3 & x80008000)>> 28;
								 case  8: r0 |= word_slice_64_16(f[j+ 7]<< 3 & x80008000)>> 32;
								 case  7: r0 |= word_slice_64_16(f[j+ 6]<< 3 & x80008000)>> 36;
								 case  6: r0 |= word_slice_64_16(f[j+ 5]<< 3 & x80008000)>> 40;
								 case  5: r0 |= word_slice_64_16(f[j+ 4]<< 3 & x80008000)>> 44;
								 case  4: r0 |= word_slice_64_16(f[j+ 3]<< 3 & x80008000)>> 48;
								 case  3: r0 |= word_slice_64_16(f[j+ 2]<< 3 & x80008000)>> 52;
								 case  2: r0 |= word_slice_64_16(f[j+ 1]<< 3 & x80008000)>> 56;
								 case  1: r0 |= word_slice_64_16(f[j+ 0]<< 3 & x80008000)>> 60;
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t0[j2] |= r0 & bitmask_end;
						 }
					 } break;
		}
	} else {
		switch(T->depth) {
			case 11: {
						 for(size_t i=0; i<T->nrows; i++) {
							 word *t0 = T->x[ 8]->rows[i];
							 word *t1 = T->x[ 9]->rows[i];
							 word *t2 = T->x[10]->rows[i];
							 const word const *f  = F->x->rows[i];

							 /* bulk of work */
							 for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
								 word_slice_64_16_combine_bulk(t0, j2, f, j,  7);
								 word_slice_64_16_combine_bulk(t1, j2, f, j,  6);
								 word_slice_64_16_combine_bulk(t2, j2, f, j,  5);
							 }
							 r0 = r1 = r2 = 0;
							 switch(F->x->width - j) {
								 case 16: r0 |= word_slice_64_16(f[j+15]<< 7 & x80008000)>>  0; r1 |= word_slice_64_16(f[j+15]<< 6 & x80008000)>>  0; r2 |= word_slice_64_16(f[j+15]<< 5 & x80008000)>>  0;
								 case 15: r0 |= word_slice_64_16(f[j+14]<< 7 & x80008000)>>  4; r1 |= word_slice_64_16(f[j+14]<< 6 & x80008000)>>  4; r2 |= word_slice_64_16(f[j+14]<< 5 & x80008000)>>  4;
								 case 14: r0 |= word_slice_64_16(f[j+13]<< 7 & x80008000)>>  8; r1 |= word_slice_64_16(f[j+13]<< 6 & x80008000)>>  8; r2 |= word_slice_64_16(f[j+13]<< 5 & x80008000)>>  8;
								 case 13: r0 |= word_slice_64_16(f[j+12]<< 7 & x80008000)>> 12; r1 |= word_slice_64_16(f[j+12]<< 6 & x80008000)>> 12; r2 |= word_slice_64_16(f[j+12]<< 5 & x80008000)>> 12;
								 case 12: r0 |= word_slice_64_16(f[j+11]<< 7 & x80008000)>> 16; r1 |= word_slice_64_16(f[j+11]<< 6 & x80008000)>> 16; r2 |= word_slice_64_16(f[j+11]<< 5 & x80008000)>> 16;
								 case 11: r0 |= word_slice_64_16(f[j+10]<< 7 & x80008000)>> 20; r1 |= word_slice_64_16(f[j+10]<< 6 & x80008000)>> 20; r2 |= word_slice_64_16(f[j+10]<< 5 & x80008000)>> 20;
								 case 10: r0 |= word_slice_64_16(f[j+ 9]<< 7 & x80008000)>> 24; r1 |= word_slice_64_16(f[j+ 9]<< 6 & x80008000)>> 24; r2 |= word_slice_64_16(f[j+ 9]<< 5 & x80008000)>> 24;
								 case  9: r0 |= word_slice_64_16(f[j+ 8]<< 7 & x80008000)>> 28; r1 |= word_slice_64_16(f[j+ 8]<< 6 & x80008000)>> 28; r2 |= word_slice_64_16(f[j+ 8]<< 5 & x80008000)>> 28;
								 case  8: r0 |= word_slice_64_16(f[j+ 7]<< 7 & x80008000)>> 32; r1 |= word_slice_64_16(f[j+ 7]<< 6 & x80008000)>> 32; r2 |= word_slice_64_16(f[j+ 7]<< 5 & x80008000)>> 32;
								 case  7: r0 |= word_slice_64_16(f[j+ 6]<< 7 & x80008000)>> 36; r1 |= word_slice_64_16(f[j+ 6]<< 6 & x80008000)>> 36; r2 |= word_slice_64_16(f[j+ 6]<< 5 & x80008000)>> 36;
								 case  6: r0 |= word_slice_64_16(f[j+ 5]<< 7 & x80008000)>> 40; r1 |= word_slice_64_16(f[j+ 5]<< 6 & x80008000)>> 40; r2 |= word_slice_64_16(f[j+ 5]<< 5 & x80008000)>> 40;
								 case  5: r0 |= word_slice_64_16(f[j+ 4]<< 7 & x80008000)>> 44; r1 |= word_slice_64_16(f[j+ 4]<< 6 & x80008000)>> 44; r2 |= word_slice_64_16(f[j+ 4]<< 5 & x80008000)>> 44;
								 case  4: r0 |= word_slice_64_16(f[j+ 3]<< 7 & x80008000)>> 48; r1 |= word_slice_64_16(f[j+ 3]<< 6 & x80008000)>> 48; r2 |= word_slice_64_16(f[j+ 3]<< 5 & x80008000)>> 48;
								 case  3: r0 |= word_slice_64_16(f[j+ 2]<< 7 & x80008000)>> 52; r1 |= word_slice_64_16(f[j+ 2]<< 6 & x80008000)>> 52; r2 |= word_slice_64_16(f[j+ 2]<< 5 & x80008000)>> 52;
								 case  2: r0 |= word_slice_64_16(f[j+ 1]<< 7 & x80008000)>> 56; r1 |= word_slice_64_16(f[j+ 1]<< 6 & x80008000)>> 56; r2 |= word_slice_64_16(f[j+ 1]<< 5 & x80008000)>> 56;
								 case  1: r0 |= word_slice_64_16(f[j+ 0]<< 7 & x80008000)>> 60; r1 |= word_slice_64_16(f[j+ 0]<< 6 & x80008000)>> 60; r2 |= word_slice_64_16(f[j+ 0]<< 5 & x80008000)>> 60;
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t0[j2] |= r0 & bitmask_end;
							 t1[j2] |= r1 & bitmask_end;
							 t2[j2] |= r2 & bitmask_end;
						 }
					 } break;
			case 10: {
						 for(size_t i=0; i<T->nrows; i++) {
							 word *t0 = T->x[ 8]->rows[i];
							 word *t1 = T->x[ 9]->rows[i];
							 const word const *f  = F->x->rows[i];

							 /* bulk of work */
							 for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
								 word_slice_64_16_combine_bulk(t0, j2, f, j,  7);
								 word_slice_64_16_combine_bulk(t1, j2, f, j,  6);
							 }
							 r0 = r1 = 0;
							 switch(F->x->width - j) {
								 case 16: r0 |= word_slice_64_16(f[j+15]<< 7 & x80008000)>>  0; r1 |= word_slice_64_16(f[j+15]<< 6 & x80008000)>>  0;
								 case 15: r0 |= word_slice_64_16(f[j+14]<< 7 & x80008000)>>  4; r1 |= word_slice_64_16(f[j+14]<< 6 & x80008000)>>  4;
								 case 14: r0 |= word_slice_64_16(f[j+13]<< 7 & x80008000)>>  8; r1 |= word_slice_64_16(f[j+13]<< 6 & x80008000)>>  8;
								 case 13: r0 |= word_slice_64_16(f[j+12]<< 7 & x80008000)>> 12; r1 |= word_slice_64_16(f[j+12]<< 6 & x80008000)>> 12;
								 case 12: r0 |= word_slice_64_16(f[j+11]<< 7 & x80008000)>> 16; r1 |= word_slice_64_16(f[j+11]<< 6 & x80008000)>> 16;
								 case 11: r0 |= word_slice_64_16(f[j+10]<< 7 & x80008000)>> 20; r1 |= word_slice_64_16(f[j+10]<< 6 & x80008000)>> 20;
								 case 10: r0 |= word_slice_64_16(f[j+ 9]<< 7 & x80008000)>> 24; r1 |= word_slice_64_16(f[j+ 9]<< 6 & x80008000)>> 24;
								 case  9: r0 |= word_slice_64_16(f[j+ 8]<< 7 & x80008000)>> 28; r1 |= word_slice_64_16(f[j+ 8]<< 6 & x80008000)>> 28;
								 case  8: r0 |= word_slice_64_16(f[j+ 7]<< 7 & x80008000)>> 32; r1 |= word_slice_64_16(f[j+ 7]<< 6 & x80008000)>> 32;
								 case  7: r0 |= word_slice_64_16(f[j+ 6]<< 7 & x80008000)>> 36; r1 |= word_slice_64_16(f[j+ 6]<< 6 & x80008000)>> 36;
								 case  6: r0 |= word_slice_64_16(f[j+ 5]<< 7 & x80008000)>> 40; r1 |= word_slice_64_16(f[j+ 5]<< 6 & x80008000)>> 40;
								 case  5: r0 |= word_slice_64_16(f[j+ 4]<< 7 & x80008000)>> 44; r1 |= word_slice_64_16(f[j+ 4]<< 6 & x80008000)>> 44;
								 case  4: r0 |= word_slice_64_16(f[j+ 3]<< 7 & x80008000)>> 48; r1 |= word_slice_64_16(f[j+ 3]<< 6 & x80008000)>> 48;
								 case  3: r0 |= word_slice_64_16(f[j+ 2]<< 7 & x80008000)>> 52; r1 |= word_slice_64_16(f[j+ 2]<< 6 & x80008000)>> 52;
								 case  2: r0 |= word_slice_64_16(f[j+ 1]<< 7 & x80008000)>> 56; r1 |= word_slice_64_16(f[j+ 1]<< 6 & x80008000)>> 56;
								 case  1: r0 |= word_slice_64_16(f[j+ 0]<< 7 & x80008000)>> 60; r1 |= word_slice_64_16(f[j+ 0]<< 6 & x80008000)>> 60;
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t0[j2] |= r0 & bitmask_end;
							 t1[j2] |= r1 & bitmask_end;
						 }
					 } break;
			case  9: {
						 for(size_t i=0; i<T->nrows; i++) {
							 word *t0 = T->x[ 8]->rows[i];
							 const word const *f  = F->x->rows[i];

							 /* bulk of work */
							 for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
								 word_slice_64_16_combine_bulk(t0, j2, f, j,  7);
							 }
							 r0 = 0;
							 switch(F->x->width - j) {
								 case 16: r0 |= word_slice_64_16(f[j+15]<< 7 & x80008000)>>  0;
								 case 15: r0 |= word_slice_64_16(f[j+14]<< 7 & x80008000)>>  4;
								 case 14: r0 |= word_slice_64_16(f[j+13]<< 7 & x80008000)>>  8;
								 case 13: r0 |= word_slice_64_16(f[j+12]<< 7 & x80008000)>> 12;
								 case 12: r0 |= word_slice_64_16(f[j+11]<< 7 & x80008000)>> 16;
								 case 11: r0 |= word_slice_64_16(f[j+10]<< 7 & x80008000)>> 20;
								 case 10: r0 |= word_slice_64_16(f[j+ 9]<< 7 & x80008000)>> 24;
								 case  9: r0 |= word_slice_64_16(f[j+ 8]<< 7 & x80008000)>> 28;
								 case  8: r0 |= word_slice_64_16(f[j+ 7]<< 7 & x80008000)>> 32;
								 case  7: r0 |= word_slice_64_16(f[j+ 6]<< 7 & x80008000)>> 36;
								 case  6: r0 |= word_slice_64_16(f[j+ 5]<< 7 & x80008000)>> 40;
								 case  5: r0 |= word_slice_64_16(f[j+ 4]<< 7 & x80008000)>> 44;
								 case  4: r0 |= word_slice_64_16(f[j+ 3]<< 7 & x80008000)>> 48;
								 case  3: r0 |= word_slice_64_16(f[j+ 2]<< 7 & x80008000)>> 52;
								 case  2: r0 |= word_slice_64_16(f[j+ 1]<< 7 & x80008000)>> 56;
								 case  1: r0 |= word_slice_64_16(f[j+ 0]<< 7 & x80008000)>> 60;
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t0[j2] |= r0 & bitmask_end;
						 }
					 } break;
			default:
					 m4ri_die("impossible");
		}
	}
	return T;
}
