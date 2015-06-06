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

static inline word word_cling_64_16(word a) {
	a = (a & xcccccccc & x__left04) | (a & xcccccccc>> 2  & x__left04)>>30;
	a = (a & xaaaaaaaa) | (a & xaaaaaaaa>> 1)>> 15;
	return a;
}

mzed_t *_mzed_cling16(mzed_t *T, const mzd_slice_t *F) {
	wi_t j,j2 = 0;

	const word bitmask_end = T->x->high_bitmask;

	if (mzd_slice_is_zero(F))
		return T;

	for(rci_t i=0; i<T->nrows; i++) {
		const word *f00 = F->x[ 0]->rows[i];
		const word *f01 = F->x[ 1]->rows[i];
		const word *f02 = F->x[ 2]->rows[i];
		const word *f03 = F->x[ 3]->rows[i];
		const word *f04 = F->x[ 4]->rows[i];
		const word *f05 = F->x[ 5]->rows[i];
		const word *f06 = F->x[ 6]->rows[i];
		const word *f07 = F->x[ 7]->rows[i];
		word *t  = T->x->rows[i];

		for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
			t[j+ 0] = (word_cling_64_16(f00[j2]<<60)>>15) | (word_cling_64_16(f01[j2]<<60)>>14) | (word_cling_64_16(f02[j2]<<60)>>13) | (word_cling_64_16(f03[j2]<<60)>>12) \
					  |       (word_cling_64_16(f04[j2]<<60)>>11) | (word_cling_64_16(f05[j2]<<60)>>10) | (word_cling_64_16(f06[j2]<<60)>> 9) | (word_cling_64_16(f07[j2]<<60)>> 8);
			t[j+ 1] = (word_cling_64_16(f00[j2]<<56)>>15) | (word_cling_64_16(f01[j2]<<56)>>14) | (word_cling_64_16(f02[j2]<<56)>>13) | (word_cling_64_16(f03[j2]<<56)>>12) \
					  |       (word_cling_64_16(f04[j2]<<56)>>11) | (word_cling_64_16(f05[j2]<<56)>>10) | (word_cling_64_16(f06[j2]<<56)>> 9) | (word_cling_64_16(f07[j2]<<56)>> 8);
			t[j+ 2] = (word_cling_64_16(f00[j2]<<52)>>15) | (word_cling_64_16(f01[j2]<<52)>>14) | (word_cling_64_16(f02[j2]<<52)>>13) | (word_cling_64_16(f03[j2]<<52)>>12) \
					  |       (word_cling_64_16(f04[j2]<<52)>>11) | (word_cling_64_16(f05[j2]<<52)>>10) | (word_cling_64_16(f06[j2]<<52)>> 9) | (word_cling_64_16(f07[j2]<<52)>> 8);
			t[j+ 3] = (word_cling_64_16(f00[j2]<<48)>>15) | (word_cling_64_16(f01[j2]<<48)>>14) | (word_cling_64_16(f02[j2]<<48)>>13) | (word_cling_64_16(f03[j2]<<48)>>12) \
					  |       (word_cling_64_16(f04[j2]<<48)>>11) | (word_cling_64_16(f05[j2]<<48)>>10) | (word_cling_64_16(f06[j2]<<48)>> 9) | (word_cling_64_16(f07[j2]<<48)>> 8);
			t[j+ 4] = (word_cling_64_16(f00[j2]<<44)>>15) | (word_cling_64_16(f01[j2]<<44)>>14) | (word_cling_64_16(f02[j2]<<44)>>13) | (word_cling_64_16(f03[j2]<<44)>>12) \
					  |       (word_cling_64_16(f04[j2]<<44)>>11) | (word_cling_64_16(f05[j2]<<44)>>10) | (word_cling_64_16(f06[j2]<<44)>> 9) | (word_cling_64_16(f07[j2]<<44)>> 8);
			t[j+ 5] = (word_cling_64_16(f00[j2]<<40)>>15) | (word_cling_64_16(f01[j2]<<40)>>14) | (word_cling_64_16(f02[j2]<<40)>>13) | (word_cling_64_16(f03[j2]<<40)>>12) \
					  |       (word_cling_64_16(f04[j2]<<40)>>11) | (word_cling_64_16(f05[j2]<<40)>>10) | (word_cling_64_16(f06[j2]<<40)>> 9) | (word_cling_64_16(f07[j2]<<40)>> 8);
			t[j+ 6] = (word_cling_64_16(f00[j2]<<36)>>15) | (word_cling_64_16(f01[j2]<<36)>>14) | (word_cling_64_16(f02[j2]<<36)>>13) | (word_cling_64_16(f03[j2]<<36)>>12) \
					  |       (word_cling_64_16(f04[j2]<<36)>>11) | (word_cling_64_16(f05[j2]<<36)>>10) | (word_cling_64_16(f06[j2]<<36)>> 9) | (word_cling_64_16(f07[j2]<<36)>> 8);
			t[j+ 7] = (word_cling_64_16(f00[j2]<<32)>>15) | (word_cling_64_16(f01[j2]<<32)>>14) | (word_cling_64_16(f02[j2]<<32)>>13) | (word_cling_64_16(f03[j2]<<32)>>12) \
					  |       (word_cling_64_16(f04[j2]<<32)>>11) | (word_cling_64_16(f05[j2]<<32)>>10) | (word_cling_64_16(f06[j2]<<32)>> 9) | (word_cling_64_16(f07[j2]<<32)>> 8);
			t[j+ 8] = (word_cling_64_16(f00[j2]<<28)>>15) | (word_cling_64_16(f01[j2]<<28)>>14) | (word_cling_64_16(f02[j2]<<28)>>13) | (word_cling_64_16(f03[j2]<<28)>>12) \
					  |       (word_cling_64_16(f04[j2]<<28)>>11) | (word_cling_64_16(f05[j2]<<28)>>10) | (word_cling_64_16(f06[j2]<<28)>> 9) | (word_cling_64_16(f07[j2]<<28)>> 8);
			t[j+ 9] = (word_cling_64_16(f00[j2]<<24)>>15) | (word_cling_64_16(f01[j2]<<24)>>14) | (word_cling_64_16(f02[j2]<<24)>>13) | (word_cling_64_16(f03[j2]<<24)>>12) \
					  |       (word_cling_64_16(f04[j2]<<24)>>11) | (word_cling_64_16(f05[j2]<<24)>>10) | (word_cling_64_16(f06[j2]<<24)>> 9) | (word_cling_64_16(f07[j2]<<24)>> 8);
			t[j+10] = (word_cling_64_16(f00[j2]<<20)>>15) | (word_cling_64_16(f01[j2]<<20)>>14) | (word_cling_64_16(f02[j2]<<20)>>13) | (word_cling_64_16(f03[j2]<<20)>>12) \
					  |       (word_cling_64_16(f04[j2]<<20)>>11) | (word_cling_64_16(f05[j2]<<20)>>10) | (word_cling_64_16(f06[j2]<<20)>> 9) | (word_cling_64_16(f07[j2]<<20)>> 8);
			t[j+11] = (word_cling_64_16(f00[j2]<<16)>>15) | (word_cling_64_16(f01[j2]<<16)>>14) | (word_cling_64_16(f02[j2]<<16)>>13) | (word_cling_64_16(f03[j2]<<16)>>12) \
					  |       (word_cling_64_16(f04[j2]<<16)>>11) | (word_cling_64_16(f05[j2]<<16)>>10) | (word_cling_64_16(f06[j2]<<16)>> 9) | (word_cling_64_16(f07[j2]<<16)>> 8);
			t[j+12] = (word_cling_64_16(f00[j2]<<12)>>15) | (word_cling_64_16(f01[j2]<<12)>>14) | (word_cling_64_16(f02[j2]<<12)>>13) | (word_cling_64_16(f03[j2]<<12)>>12) \
					  |       (word_cling_64_16(f04[j2]<<12)>>11) | (word_cling_64_16(f05[j2]<<12)>>10) | (word_cling_64_16(f06[j2]<<12)>> 9) | (word_cling_64_16(f07[j2]<<12)>> 8);
			t[j+13] = (word_cling_64_16(f00[j2]<< 8)>>15) | (word_cling_64_16(f01[j2]<< 8)>>14) | (word_cling_64_16(f02[j2]<< 8)>>13) | (word_cling_64_16(f03[j2]<< 8)>>12) \
					  |       (word_cling_64_16(f04[j2]<< 8)>>11) | (word_cling_64_16(f05[j2]<< 8)>>10) | (word_cling_64_16(f06[j2]<< 8)>> 9) | (word_cling_64_16(f07[j2]<< 8)>> 8);
			t[j+14] = (word_cling_64_16(f00[j2]<< 4)>>15) | (word_cling_64_16(f01[j2]<< 4)>>14) | (word_cling_64_16(f02[j2]<< 4)>>13) | (word_cling_64_16(f03[j2]<< 4)>>12) \
					  |       (word_cling_64_16(f04[j2]<< 4)>>11) | (word_cling_64_16(f05[j2]<< 4)>>10) | (word_cling_64_16(f06[j2]<< 4)>> 9) | (word_cling_64_16(f07[j2]<< 4)>> 8);
			t[j+15] = (word_cling_64_16(f00[j2]<< 0)>>15) | (word_cling_64_16(f01[j2]<< 0)>>14) | (word_cling_64_16(f02[j2]<< 0)>>13) | (word_cling_64_16(f03[j2]<< 0)>>12) \
					  |       (word_cling_64_16(f04[j2]<< 0)>>11) | (word_cling_64_16(f05[j2]<< 0)>>10) | (word_cling_64_16(f06[j2]<< 0)>> 9) | (word_cling_64_16(f07[j2]<< 0)>> 8);
		}

		register word tmp = t[T->x->width-1];
		switch(T->x->width - j) {
			case 16: t[j+15] = (word_cling_64_16(f00[j2]<< 0)>>15) | (word_cling_64_16(f01[j2]<< 0)>>14) | (word_cling_64_16(f02[j2]<< 0)>>13) | (word_cling_64_16(f03[j2]<< 0)>>12) | \
					 (word_cling_64_16(f04[j2]<< 0)>>11) | (word_cling_64_16(f05[j2]<< 0)>>10) | (word_cling_64_16(f06[j2]<< 0)>> 9) | (word_cling_64_16(f07[j2]<< 0)>> 8);
			case 15: t[j+14] = (word_cling_64_16(f00[j2]<< 4)>>15) | (word_cling_64_16(f01[j2]<< 4)>>14) | (word_cling_64_16(f02[j2]<< 4)>>13) | (word_cling_64_16(f03[j2]<< 4)>>12) | \
					 (word_cling_64_16(f04[j2]<< 4)>>11) | (word_cling_64_16(f05[j2]<< 4)>>10) | (word_cling_64_16(f06[j2]<< 4)>> 9) | (word_cling_64_16(f07[j2]<< 4)>> 8);
			case 14: t[j+13] = (word_cling_64_16(f00[j2]<< 8)>>15) | (word_cling_64_16(f01[j2]<< 8)>>14) | (word_cling_64_16(f02[j2]<< 8)>>13) | (word_cling_64_16(f03[j2]<< 8)>>12) | \
					 (word_cling_64_16(f04[j2]<< 8)>>11) | (word_cling_64_16(f05[j2]<< 8)>>10) | (word_cling_64_16(f06[j2]<< 8)>> 9) | (word_cling_64_16(f07[j2]<< 8)>> 8);
			case 13: t[j+12] = (word_cling_64_16(f00[j2]<<12)>>15) | (word_cling_64_16(f01[j2]<<12)>>14) | (word_cling_64_16(f02[j2]<<12)>>13) | (word_cling_64_16(f03[j2]<<12)>>12) | \
					 (word_cling_64_16(f04[j2]<<12)>>11) | (word_cling_64_16(f05[j2]<<12)>>10) | (word_cling_64_16(f06[j2]<<12)>> 9) | (word_cling_64_16(f07[j2]<<12)>> 8);
			case 12: t[j+11] = (word_cling_64_16(f00[j2]<<16)>>15) | (word_cling_64_16(f01[j2]<<16)>>14) | (word_cling_64_16(f02[j2]<<16)>>13) | (word_cling_64_16(f03[j2]<<16)>>12) | \
					 (word_cling_64_16(f04[j2]<<16)>>11) | (word_cling_64_16(f05[j2]<<16)>>10) | (word_cling_64_16(f06[j2]<<16)>> 9) | (word_cling_64_16(f07[j2]<<16)>> 8);
			case 11: t[j+10] = (word_cling_64_16(f00[j2]<<20)>>15) | (word_cling_64_16(f01[j2]<<20)>>14) | (word_cling_64_16(f02[j2]<<20)>>13) | (word_cling_64_16(f03[j2]<<20)>>12) | \
					 (word_cling_64_16(f04[j2]<<20)>>11) | (word_cling_64_16(f05[j2]<<20)>>10) | (word_cling_64_16(f06[j2]<<20)>> 9) | (word_cling_64_16(f07[j2]<<20)>> 8);
			case 10: t[j+ 9] = (word_cling_64_16(f00[j2]<<24)>>15) | (word_cling_64_16(f01[j2]<<24)>>14) | (word_cling_64_16(f02[j2]<<24)>>13) | (word_cling_64_16(f03[j2]<<24)>>12) | \
					 (word_cling_64_16(f04[j2]<<24)>>11) | (word_cling_64_16(f05[j2]<<24)>>10) | (word_cling_64_16(f06[j2]<<24)>> 9) | (word_cling_64_16(f07[j2]<<24)>> 8);
			case  9: t[j+ 8] = (word_cling_64_16(f00[j2]<<28)>>15) | (word_cling_64_16(f01[j2]<<28)>>14) | (word_cling_64_16(f02[j2]<<28)>>13) | (word_cling_64_16(f03[j2]<<28)>>12) | \
					 (word_cling_64_16(f04[j2]<<28)>>11) | (word_cling_64_16(f05[j2]<<28)>>10) | (word_cling_64_16(f06[j2]<<28)>> 9) | (word_cling_64_16(f07[j2]<<28)>> 8);
			case  8: t[j+ 7] = (word_cling_64_16(f00[j2]<<32)>>15) | (word_cling_64_16(f01[j2]<<32)>>14) | (word_cling_64_16(f02[j2]<<32)>>13) | (word_cling_64_16(f03[j2]<<32)>>12) | \
					 (word_cling_64_16(f04[j2]<<32)>>11) | (word_cling_64_16(f05[j2]<<32)>>10) | (word_cling_64_16(f06[j2]<<32)>> 9) | (word_cling_64_16(f07[j2]<<32)>> 8);
			case  7: t[j+ 6] = (word_cling_64_16(f00[j2]<<36)>>15) | (word_cling_64_16(f01[j2]<<36)>>14) | (word_cling_64_16(f02[j2]<<36)>>13) | (word_cling_64_16(f03[j2]<<36)>>12) | \
					 (word_cling_64_16(f04[j2]<<36)>>11) | (word_cling_64_16(f05[j2]<<36)>>10) | (word_cling_64_16(f06[j2]<<36)>> 9) | (word_cling_64_16(f07[j2]<<36)>> 8);
			case  6: t[j+ 5] = (word_cling_64_16(f00[j2]<<40)>>15) | (word_cling_64_16(f01[j2]<<40)>>14) | (word_cling_64_16(f02[j2]<<40)>>13) | (word_cling_64_16(f03[j2]<<40)>>12) | \
					 (word_cling_64_16(f04[j2]<<40)>>11) | (word_cling_64_16(f05[j2]<<40)>>10) | (word_cling_64_16(f06[j2]<<40)>> 9) | (word_cling_64_16(f07[j2]<<40)>> 8);
			case  5: t[j+ 4] = (word_cling_64_16(f00[j2]<<44)>>15) | (word_cling_64_16(f01[j2]<<44)>>14) | (word_cling_64_16(f02[j2]<<44)>>13) | (word_cling_64_16(f03[j2]<<44)>>12) | \
					 (word_cling_64_16(f04[j2]<<44)>>11) | (word_cling_64_16(f05[j2]<<44)>>10) | (word_cling_64_16(f06[j2]<<44)>> 9) | (word_cling_64_16(f07[j2]<<44)>> 8);
			case  4: t[j+ 3] = (word_cling_64_16(f00[j2]<<48)>>15) | (word_cling_64_16(f01[j2]<<48)>>14) | (word_cling_64_16(f02[j2]<<48)>>13) | (word_cling_64_16(f03[j2]<<48)>>12) | \
					 (word_cling_64_16(f04[j2]<<48)>>11) | (word_cling_64_16(f05[j2]<<48)>>10) | (word_cling_64_16(f06[j2]<<48)>> 9) | (word_cling_64_16(f07[j2]<<48)>> 8);
			case  3: t[j+ 2] = (word_cling_64_16(f00[j2]<<52)>>15) | (word_cling_64_16(f01[j2]<<52)>>14) | (word_cling_64_16(f02[j2]<<52)>>13) | (word_cling_64_16(f03[j2]<<52)>>12) | \
					 (word_cling_64_16(f04[j2]<<52)>>11) | (word_cling_64_16(f05[j2]<<52)>>10) | (word_cling_64_16(f06[j2]<<52)>> 9) | (word_cling_64_16(f07[j2]<<52)>> 8);
			case  2: t[j+ 1] = (word_cling_64_16(f00[j2]<<56)>>15) | (word_cling_64_16(f01[j2]<<56)>>14) | (word_cling_64_16(f02[j2]<<56)>>13) | (word_cling_64_16(f03[j2]<<56)>>12) | \
					 (word_cling_64_16(f04[j2]<<56)>>11) | (word_cling_64_16(f05[j2]<<56)>>10) | (word_cling_64_16(f06[j2]<<56)>> 9) | (word_cling_64_16(f07[j2]<<56)>> 8);
			case  1: t[j+ 0] = (word_cling_64_16(f00[j2]<<60)>>15) | (word_cling_64_16(f01[j2]<<60)>>14) | (word_cling_64_16(f02[j2]<<60)>>13) | (word_cling_64_16(f03[j2]<<60)>>12) | \
					 (word_cling_64_16(f04[j2]<<60)>>11) | (word_cling_64_16(f05[j2]<<60)>>10) | (word_cling_64_16(f06[j2]<<60)>> 9) | (word_cling_64_16(f07[j2]<<60)>> 8);
					 break;
			default:
					 m4ri_die("impossible");
		}
		t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
	}

	if(T->finite_field->degree < 12) {
		switch(T->finite_field->degree) {
			case 9: {
						for(rci_t i=0; i<T->nrows; i++) {
							const word *f00 = F->x[ 8]->rows[i];
							word *t  = T->x->rows[i];

							for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
								t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>7);
								t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>7);
								t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>7);
								t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>7);
								t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>7);
								t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>7);
								t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>7);
								t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>7);
								t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>7);
								t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>7);
								t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>7);
								t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>7);
								t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>7);
								t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>7);
								t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>7);
								t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>7);
							}

							register word tmp = t[T->x->width-1];
							switch(T->x->width - j) {
								case 16: t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>7);
								case 15: t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>7);
								case 14: t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>7);
								case 13: t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>7);
								case 12: t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>7);
								case 11: t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>7);
								case 10: t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>7);
								case  9: t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>7);
								case  8: t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>7);
								case  7: t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>7);
								case  6: t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>7);
								case  5: t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>7);
								case  4: t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>7);
								case  3: t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>7);
								case  2: t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>7);
								case  1: t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>7);
										 break;
								default:
										 m4ri_die("impossible");
							}
							t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
						}
					}
					break;
			case 10: {
						 for(rci_t i=0; i<T->nrows; i++) {
							 const word *f00 = F->x[ 8]->rows[i];
							 const word *f01 = F->x[ 9]->rows[i];
							 word *t  = T->x->rows[i];

							 for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
								 t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>7) | (word_cling_64_16(f01[j2]<<60)>>6);
								 t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>7) | (word_cling_64_16(f01[j2]<<56)>>6);
								 t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>7) | (word_cling_64_16(f01[j2]<<52)>>6);
								 t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>7) | (word_cling_64_16(f01[j2]<<48)>>6);
								 t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>7) | (word_cling_64_16(f01[j2]<<44)>>6);
								 t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>7) | (word_cling_64_16(f01[j2]<<40)>>6);
								 t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>7) | (word_cling_64_16(f01[j2]<<36)>>6);
								 t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>7) | (word_cling_64_16(f01[j2]<<32)>>6);
								 t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>7) | (word_cling_64_16(f01[j2]<<28)>>6);
								 t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>7) | (word_cling_64_16(f01[j2]<<24)>>6);
								 t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>7) | (word_cling_64_16(f01[j2]<<20)>>6);
								 t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>7) | (word_cling_64_16(f01[j2]<<16)>>6);
								 t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>7) | (word_cling_64_16(f01[j2]<<12)>>6);
								 t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>7) | (word_cling_64_16(f01[j2]<< 8)>>6);
								 t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>7) | (word_cling_64_16(f01[j2]<< 4)>>6);
								 t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>7) | (word_cling_64_16(f01[j2]<< 0)>>6);
							 }

							 register word tmp = t[T->x->width-1];
							 switch(T->x->width - j) {
								 case 16: t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>7) | (word_cling_64_16(f01[j2]<< 0)>>6);
								 case 15: t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>7) | (word_cling_64_16(f01[j2]<< 4)>>6);
								 case 14: t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>7) | (word_cling_64_16(f01[j2]<< 8)>>6);
								 case 13: t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>7) | (word_cling_64_16(f01[j2]<<12)>>6);
								 case 12: t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>7) | (word_cling_64_16(f01[j2]<<16)>>6);
								 case 11: t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>7) | (word_cling_64_16(f01[j2]<<20)>>6);
								 case 10: t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>7) | (word_cling_64_16(f01[j2]<<24)>>6);
								 case  9: t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>7) | (word_cling_64_16(f01[j2]<<28)>>6);
								 case  8: t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>7) | (word_cling_64_16(f01[j2]<<32)>>6);
								 case  7: t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>7) | (word_cling_64_16(f01[j2]<<36)>>6);
								 case  6: t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>7) | (word_cling_64_16(f01[j2]<<40)>>6);
								 case  5: t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>7) | (word_cling_64_16(f01[j2]<<44)>>6);
								 case  4: t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>7) | (word_cling_64_16(f01[j2]<<48)>>6);
								 case  3: t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>7) | (word_cling_64_16(f01[j2]<<52)>>6);
								 case  2: t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>7) | (word_cling_64_16(f01[j2]<<56)>>6);
								 case  1: t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>7) | (word_cling_64_16(f01[j2]<<60)>>6);
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
						 }
					 }
					 break;
			case 11: {
						 for(rci_t i=0; i<T->nrows; i++) {
							 const word *f00 = F->x[ 8]->rows[i];
							 const word *f01 = F->x[ 9]->rows[i];
							 const word *f02 = F->x[10]->rows[i];
							 word *t  = T->x->rows[i];

							 for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
								 t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>7) | (word_cling_64_16(f01[j2]<<60)>>6) | (word_cling_64_16(f02[j2]<<60)>>5);
								 t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>7) | (word_cling_64_16(f01[j2]<<56)>>6) | (word_cling_64_16(f02[j2]<<56)>>5);
								 t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>7) | (word_cling_64_16(f01[j2]<<52)>>6) | (word_cling_64_16(f02[j2]<<52)>>5);
								 t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>7) | (word_cling_64_16(f01[j2]<<48)>>6) | (word_cling_64_16(f02[j2]<<48)>>5);
								 t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>7) | (word_cling_64_16(f01[j2]<<44)>>6) | (word_cling_64_16(f02[j2]<<44)>>5);
								 t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>7) | (word_cling_64_16(f01[j2]<<40)>>6) | (word_cling_64_16(f02[j2]<<40)>>5);
								 t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>7) | (word_cling_64_16(f01[j2]<<36)>>6) | (word_cling_64_16(f02[j2]<<36)>>5);
								 t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>7) | (word_cling_64_16(f01[j2]<<32)>>6) | (word_cling_64_16(f02[j2]<<32)>>5);
								 t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>7) | (word_cling_64_16(f01[j2]<<28)>>6) | (word_cling_64_16(f02[j2]<<28)>>5);
								 t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>7) | (word_cling_64_16(f01[j2]<<24)>>6) | (word_cling_64_16(f02[j2]<<24)>>5);
								 t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>7) | (word_cling_64_16(f01[j2]<<20)>>6) | (word_cling_64_16(f02[j2]<<20)>>5);
								 t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>7) | (word_cling_64_16(f01[j2]<<16)>>6) | (word_cling_64_16(f02[j2]<<16)>>5);
								 t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>7) | (word_cling_64_16(f01[j2]<<12)>>6) | (word_cling_64_16(f02[j2]<<12)>>5);
								 t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>7) | (word_cling_64_16(f01[j2]<< 8)>>6) | (word_cling_64_16(f02[j2]<< 8)>>5);
								 t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>7) | (word_cling_64_16(f01[j2]<< 4)>>6) | (word_cling_64_16(f02[j2]<< 4)>>5);
								 t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>7) | (word_cling_64_16(f01[j2]<< 0)>>6) | (word_cling_64_16(f02[j2]<< 0)>>5);
							 }

							 register word tmp = t[T->x->width-1];
							 switch(T->x->width - j) {
								 case 16: t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>7) | (word_cling_64_16(f01[j2]<< 0)>>6) | (word_cling_64_16(f02[j2]<< 0)>>5);
								 case 15: t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>7) | (word_cling_64_16(f01[j2]<< 4)>>6) | (word_cling_64_16(f02[j2]<< 4)>>5);
								 case 14: t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>7) | (word_cling_64_16(f01[j2]<< 8)>>6) | (word_cling_64_16(f02[j2]<< 8)>>5);
								 case 13: t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>7) | (word_cling_64_16(f01[j2]<<12)>>6) | (word_cling_64_16(f02[j2]<<12)>>5);
								 case 12: t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>7) | (word_cling_64_16(f01[j2]<<16)>>6) | (word_cling_64_16(f02[j2]<<16)>>5);
								 case 11: t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>7) | (word_cling_64_16(f01[j2]<<20)>>6) | (word_cling_64_16(f02[j2]<<20)>>5);
								 case 10: t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>7) | (word_cling_64_16(f01[j2]<<24)>>6) | (word_cling_64_16(f02[j2]<<24)>>5);
								 case  9: t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>7) | (word_cling_64_16(f01[j2]<<28)>>6) | (word_cling_64_16(f02[j2]<<28)>>5);
								 case  8: t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>7) | (word_cling_64_16(f01[j2]<<32)>>6) | (word_cling_64_16(f02[j2]<<32)>>5);
								 case  7: t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>7) | (word_cling_64_16(f01[j2]<<36)>>6) | (word_cling_64_16(f02[j2]<<36)>>5);
								 case  6: t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>7) | (word_cling_64_16(f01[j2]<<40)>>6) | (word_cling_64_16(f02[j2]<<40)>>5);
								 case  5: t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>7) | (word_cling_64_16(f01[j2]<<44)>>6) | (word_cling_64_16(f02[j2]<<44)>>5);
								 case  4: t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>7) | (word_cling_64_16(f01[j2]<<48)>>6) | (word_cling_64_16(f02[j2]<<48)>>5);
								 case  3: t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>7) | (word_cling_64_16(f01[j2]<<52)>>6) | (word_cling_64_16(f02[j2]<<52)>>5);
								 case  2: t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>7) | (word_cling_64_16(f01[j2]<<56)>>6) | (word_cling_64_16(f02[j2]<<56)>>5);
								 case  1: t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>7) | (word_cling_64_16(f01[j2]<<60)>>6) | (word_cling_64_16(f02[j2]<<60)>>5);
										  break;
								 default:
										  m4ri_die("impossible");
							 }
							 t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
						 }
					 }
					 break;
		}
	} else {
		for(rci_t i=0; i<T->nrows; i++) {
			const word *f00 = F->x[ 8]->rows[i];
			const word *f01 = F->x[ 9]->rows[i];
			const word *f02 = F->x[10]->rows[i];
			const word *f03 = F->x[11]->rows[i];
			word *t  = T->x->rows[i];

			for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
				t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>7) | (word_cling_64_16(f01[j2]<<60)>>6) | (word_cling_64_16(f02[j2]<<60)>>5) | (word_cling_64_16(f03[j2]<<60)>>4);
				t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>7) | (word_cling_64_16(f01[j2]<<56)>>6) | (word_cling_64_16(f02[j2]<<56)>>5) | (word_cling_64_16(f03[j2]<<56)>>4);
				t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>7) | (word_cling_64_16(f01[j2]<<52)>>6) | (word_cling_64_16(f02[j2]<<52)>>5) | (word_cling_64_16(f03[j2]<<52)>>4);
				t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>7) | (word_cling_64_16(f01[j2]<<48)>>6) | (word_cling_64_16(f02[j2]<<48)>>5) | (word_cling_64_16(f03[j2]<<48)>>4);
				t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>7) | (word_cling_64_16(f01[j2]<<44)>>6) | (word_cling_64_16(f02[j2]<<44)>>5) | (word_cling_64_16(f03[j2]<<44)>>4);
				t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>7) | (word_cling_64_16(f01[j2]<<40)>>6) | (word_cling_64_16(f02[j2]<<40)>>5) | (word_cling_64_16(f03[j2]<<40)>>4);
				t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>7) | (word_cling_64_16(f01[j2]<<36)>>6) | (word_cling_64_16(f02[j2]<<36)>>5) | (word_cling_64_16(f03[j2]<<36)>>4);
				t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>7) | (word_cling_64_16(f01[j2]<<32)>>6) | (word_cling_64_16(f02[j2]<<32)>>5) | (word_cling_64_16(f03[j2]<<32)>>4);
				t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>7) | (word_cling_64_16(f01[j2]<<28)>>6) | (word_cling_64_16(f02[j2]<<28)>>5) | (word_cling_64_16(f03[j2]<<28)>>4);
				t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>7) | (word_cling_64_16(f01[j2]<<24)>>6) | (word_cling_64_16(f02[j2]<<24)>>5) | (word_cling_64_16(f03[j2]<<24)>>4);
				t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>7) | (word_cling_64_16(f01[j2]<<20)>>6) | (word_cling_64_16(f02[j2]<<20)>>5) | (word_cling_64_16(f03[j2]<<20)>>4);
				t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>7) | (word_cling_64_16(f01[j2]<<16)>>6) | (word_cling_64_16(f02[j2]<<16)>>5) | (word_cling_64_16(f03[j2]<<16)>>4);
				t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>7) | (word_cling_64_16(f01[j2]<<12)>>6) | (word_cling_64_16(f02[j2]<<12)>>5) | (word_cling_64_16(f03[j2]<<12)>>4);
				t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>7) | (word_cling_64_16(f01[j2]<< 8)>>6) | (word_cling_64_16(f02[j2]<< 8)>>5) | (word_cling_64_16(f03[j2]<< 8)>>4);
				t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>7) | (word_cling_64_16(f01[j2]<< 4)>>6) | (word_cling_64_16(f02[j2]<< 4)>>5) | (word_cling_64_16(f03[j2]<< 4)>>4);
				t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>7) | (word_cling_64_16(f01[j2]<< 0)>>6) | (word_cling_64_16(f02[j2]<< 0)>>5) | (word_cling_64_16(f03[j2]<< 0)>>4);
			}

			register word tmp = t[T->x->width-1];
			switch(T->x->width - j) {
				case 16: t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>7) | (word_cling_64_16(f01[j2]<< 0)>>6) | (word_cling_64_16(f02[j2]<< 0)>>5) | (word_cling_64_16(f03[j2]<< 0)>>4);
				case 15: t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>7) | (word_cling_64_16(f01[j2]<< 4)>>6) | (word_cling_64_16(f02[j2]<< 4)>>5) | (word_cling_64_16(f03[j2]<< 4)>>4);
				case 14: t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>7) | (word_cling_64_16(f01[j2]<< 8)>>6) | (word_cling_64_16(f02[j2]<< 8)>>5) | (word_cling_64_16(f03[j2]<< 8)>>4);
				case 13: t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>7) | (word_cling_64_16(f01[j2]<<12)>>6) | (word_cling_64_16(f02[j2]<<12)>>5) | (word_cling_64_16(f03[j2]<<12)>>4);
				case 12: t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>7) | (word_cling_64_16(f01[j2]<<16)>>6) | (word_cling_64_16(f02[j2]<<16)>>5) | (word_cling_64_16(f03[j2]<<16)>>4);
				case 11: t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>7) | (word_cling_64_16(f01[j2]<<20)>>6) | (word_cling_64_16(f02[j2]<<20)>>5) | (word_cling_64_16(f03[j2]<<20)>>4);
				case 10: t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>7) | (word_cling_64_16(f01[j2]<<24)>>6) | (word_cling_64_16(f02[j2]<<24)>>5) | (word_cling_64_16(f03[j2]<<24)>>4);
				case  9: t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>7) | (word_cling_64_16(f01[j2]<<28)>>6) | (word_cling_64_16(f02[j2]<<28)>>5) | (word_cling_64_16(f03[j2]<<28)>>4);
				case  8: t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>7) | (word_cling_64_16(f01[j2]<<32)>>6) | (word_cling_64_16(f02[j2]<<32)>>5) | (word_cling_64_16(f03[j2]<<32)>>4);
				case  7: t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>7) | (word_cling_64_16(f01[j2]<<36)>>6) | (word_cling_64_16(f02[j2]<<36)>>5) | (word_cling_64_16(f03[j2]<<36)>>4);
				case  6: t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>7) | (word_cling_64_16(f01[j2]<<40)>>6) | (word_cling_64_16(f02[j2]<<40)>>5) | (word_cling_64_16(f03[j2]<<40)>>4);
				case  5: t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>7) | (word_cling_64_16(f01[j2]<<44)>>6) | (word_cling_64_16(f02[j2]<<44)>>5) | (word_cling_64_16(f03[j2]<<44)>>4);
				case  4: t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>7) | (word_cling_64_16(f01[j2]<<48)>>6) | (word_cling_64_16(f02[j2]<<48)>>5) | (word_cling_64_16(f03[j2]<<48)>>4);
				case  3: t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>7) | (word_cling_64_16(f01[j2]<<52)>>6) | (word_cling_64_16(f02[j2]<<52)>>5) | (word_cling_64_16(f03[j2]<<52)>>4);
				case  2: t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>7) | (word_cling_64_16(f01[j2]<<56)>>6) | (word_cling_64_16(f02[j2]<<56)>>5) | (word_cling_64_16(f03[j2]<<56)>>4);
				case  1: t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>7) | (word_cling_64_16(f01[j2]<<60)>>6) | (word_cling_64_16(f02[j2]<<60)>>5) | (word_cling_64_16(f03[j2]<<60)>>4);
						 break;
				default:
						 m4ri_die("impossible");
			}
			t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);

			switch(T->finite_field->degree) {
				case 13: {
							 for(rci_t i=0; i<T->nrows; i++) {
								 const word *f00 = F->x[12]->rows[i];
								 word *t  = T->x->rows[i];

								 for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
									 t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>3);
									 t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>3);
									 t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>3);
									 t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>3);
									 t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>3);
									 t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>3);
									 t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>3);
									 t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>3);
									 t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>3);
									 t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>3);
									 t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>3);
									 t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>3);
									 t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>3);
									 t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>3);
									 t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>3);
									 t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>3);
								 }

								 register word tmp = t[T->x->width-1];
								 switch(T->x->width - j) {
									 case 16: t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>3);
									 case 15: t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>3);
									 case 14: t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>3);
									 case 13: t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>3);
									 case 12: t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>3);
									 case 11: t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>3);
									 case 10: t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>3);
									 case  9: t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>3);
									 case  8: t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>3);
									 case  7: t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>3);
									 case  6: t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>3);
									 case  5: t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>3);
									 case  4: t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>3);
									 case  3: t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>3);
									 case  2: t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>3);
									 case  1: t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>3);
											  break;
									 default:
											  m4ri_die("impossible");
								 }
								 t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
							 }
						 }
						 break;
				case 14: {
							 for(rci_t i=0; i<T->nrows; i++) {
								 const word *f00 = F->x[12]->rows[i];
								 const word *f01 = F->x[13]->rows[i];
								 word *t  = T->x->rows[i];

								 for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
									 t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>3) | (word_cling_64_16(f01[j2]<<60)>>2);
									 t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>3) | (word_cling_64_16(f01[j2]<<56)>>2);
									 t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>3) | (word_cling_64_16(f01[j2]<<52)>>2);
									 t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>3) | (word_cling_64_16(f01[j2]<<48)>>2);
									 t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>3) | (word_cling_64_16(f01[j2]<<44)>>2);
									 t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>3) | (word_cling_64_16(f01[j2]<<40)>>2);
									 t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>3) | (word_cling_64_16(f01[j2]<<36)>>2);
									 t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>3) | (word_cling_64_16(f01[j2]<<32)>>2);
									 t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>3) | (word_cling_64_16(f01[j2]<<28)>>2);
									 t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>3) | (word_cling_64_16(f01[j2]<<24)>>2);
									 t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>3) | (word_cling_64_16(f01[j2]<<20)>>2);
									 t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>3) | (word_cling_64_16(f01[j2]<<16)>>2);
									 t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>3) | (word_cling_64_16(f01[j2]<<12)>>2);
									 t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>3) | (word_cling_64_16(f01[j2]<< 8)>>2);
									 t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>3) | (word_cling_64_16(f01[j2]<< 4)>>2);
									 t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>3) | (word_cling_64_16(f01[j2]<< 0)>>2);
								 }

								 register word tmp = t[T->x->width-1];
								 switch(T->x->width - j) {
									 case 16: t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>3) | (word_cling_64_16(f01[j2]<< 0)>>2);
									 case 15: t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>3) | (word_cling_64_16(f01[j2]<< 4)>>2);
									 case 14: t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>3) | (word_cling_64_16(f01[j2]<< 8)>>2);
									 case 13: t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>3) | (word_cling_64_16(f01[j2]<<12)>>2);
									 case 12: t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>3) | (word_cling_64_16(f01[j2]<<16)>>2);
									 case 11: t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>3) | (word_cling_64_16(f01[j2]<<20)>>2);
									 case 10: t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>3) | (word_cling_64_16(f01[j2]<<24)>>2);
									 case  9: t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>3) | (word_cling_64_16(f01[j2]<<28)>>2);
									 case  8: t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>3) | (word_cling_64_16(f01[j2]<<32)>>2);
									 case  7: t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>3) | (word_cling_64_16(f01[j2]<<36)>>2);
									 case  6: t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>3) | (word_cling_64_16(f01[j2]<<40)>>2);
									 case  5: t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>3) | (word_cling_64_16(f01[j2]<<44)>>2);
									 case  4: t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>3) | (word_cling_64_16(f01[j2]<<48)>>2);
									 case  3: t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>3) | (word_cling_64_16(f01[j2]<<52)>>2);
									 case  2: t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>3) | (word_cling_64_16(f01[j2]<<56)>>2);
									 case  1: t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>3) | (word_cling_64_16(f01[j2]<<60)>>2);
											  break;
									 default:
											  m4ri_die("impossible");
								 }
								 t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
							 }
						 }
						 break;
				case 15: {
							 for(rci_t i=0; i<T->nrows; i++) {
								 const word *f00 = F->x[12]->rows[i];
								 const word *f01 = F->x[13]->rows[i];
								 const word *f02 = F->x[14]->rows[i];
								 word *t  = T->x->rows[i];

								 for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
									 t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>3) | (word_cling_64_16(f01[j2]<<60)>>2) | (word_cling_64_16(f02[j2]<<60)>>1);
									 t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>3) | (word_cling_64_16(f01[j2]<<56)>>2) | (word_cling_64_16(f02[j2]<<56)>>1);
									 t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>3) | (word_cling_64_16(f01[j2]<<52)>>2) | (word_cling_64_16(f02[j2]<<52)>>1);
									 t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>3) | (word_cling_64_16(f01[j2]<<48)>>2) | (word_cling_64_16(f02[j2]<<48)>>1);
									 t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>3) | (word_cling_64_16(f01[j2]<<44)>>2) | (word_cling_64_16(f02[j2]<<44)>>1);
									 t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>3) | (word_cling_64_16(f01[j2]<<40)>>2) | (word_cling_64_16(f02[j2]<<40)>>1);
									 t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>3) | (word_cling_64_16(f01[j2]<<36)>>2) | (word_cling_64_16(f02[j2]<<36)>>1);
									 t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>3) | (word_cling_64_16(f01[j2]<<32)>>2) | (word_cling_64_16(f02[j2]<<32)>>1);
									 t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>3) | (word_cling_64_16(f01[j2]<<28)>>2) | (word_cling_64_16(f02[j2]<<28)>>1);
									 t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>3) | (word_cling_64_16(f01[j2]<<24)>>2) | (word_cling_64_16(f02[j2]<<24)>>1);
									 t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>3) | (word_cling_64_16(f01[j2]<<20)>>2) | (word_cling_64_16(f02[j2]<<20)>>1);
									 t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>3) | (word_cling_64_16(f01[j2]<<16)>>2) | (word_cling_64_16(f02[j2]<<16)>>1);
									 t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>3) | (word_cling_64_16(f01[j2]<<12)>>2) | (word_cling_64_16(f02[j2]<<12)>>1);
									 t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>3) | (word_cling_64_16(f01[j2]<< 8)>>2) | (word_cling_64_16(f02[j2]<< 8)>>1);
									 t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>3) | (word_cling_64_16(f01[j2]<< 4)>>2) | (word_cling_64_16(f02[j2]<< 4)>>1);
									 t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>3) | (word_cling_64_16(f01[j2]<< 0)>>2) | (word_cling_64_16(f02[j2]<< 0)>>1);
								 }

								 register word tmp = t[T->x->width-1];
								 switch(T->x->width - j) {
									 case 16: t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>3) | (word_cling_64_16(f01[j2]<< 0)>>2) | (word_cling_64_16(f02[j2]<< 0)>>1);
									 case 15: t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>3) | (word_cling_64_16(f01[j2]<< 4)>>2) | (word_cling_64_16(f02[j2]<< 4)>>1);
									 case 14: t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>3) | (word_cling_64_16(f01[j2]<< 8)>>2) | (word_cling_64_16(f02[j2]<< 8)>>1);
									 case 13: t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>3) | (word_cling_64_16(f01[j2]<<12)>>2) | (word_cling_64_16(f02[j2]<<12)>>1);
									 case 12: t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>3) | (word_cling_64_16(f01[j2]<<16)>>2) | (word_cling_64_16(f02[j2]<<16)>>1);
									 case 11: t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>3) | (word_cling_64_16(f01[j2]<<20)>>2) | (word_cling_64_16(f02[j2]<<20)>>1);
									 case 10: t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>3) | (word_cling_64_16(f01[j2]<<24)>>2) | (word_cling_64_16(f02[j2]<<24)>>1);
									 case  9: t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>3) | (word_cling_64_16(f01[j2]<<28)>>2) | (word_cling_64_16(f02[j2]<<28)>>1);
									 case  8: t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>3) | (word_cling_64_16(f01[j2]<<32)>>2) | (word_cling_64_16(f02[j2]<<32)>>1);
									 case  7: t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>3) | (word_cling_64_16(f01[j2]<<36)>>2) | (word_cling_64_16(f02[j2]<<36)>>1);
									 case  6: t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>3) | (word_cling_64_16(f01[j2]<<40)>>2) | (word_cling_64_16(f02[j2]<<40)>>1);
									 case  5: t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>3) | (word_cling_64_16(f01[j2]<<44)>>2) | (word_cling_64_16(f02[j2]<<44)>>1);
									 case  4: t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>3) | (word_cling_64_16(f01[j2]<<48)>>2) | (word_cling_64_16(f02[j2]<<48)>>1);
									 case  3: t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>3) | (word_cling_64_16(f01[j2]<<52)>>2) | (word_cling_64_16(f02[j2]<<52)>>1);
									 case  2: t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>3) | (word_cling_64_16(f01[j2]<<56)>>2) | (word_cling_64_16(f02[j2]<<56)>>1);
									 case  1: t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>3) | (word_cling_64_16(f01[j2]<<60)>>2) | (word_cling_64_16(f02[j2]<<60)>>1);
											  break;
									 default:
											  m4ri_die("impossible");
								 }
								 t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
							 }
						 }
						 break;
				case 16: {
							 for(rci_t i=0; i<T->nrows; i++) {
								 const word *f00 = F->x[12]->rows[i];
								 const word *f01 = F->x[13]->rows[i];
								 const word *f02 = F->x[14]->rows[i];
								 const word *f03 = F->x[15]->rows[i];
								 word *t  = T->x->rows[i];

								 for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
									 t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>3) | (word_cling_64_16(f01[j2]<<60)>>2) | (word_cling_64_16(f02[j2]<<60)>>1) | (word_cling_64_16(f03[j2]<<60)>>0);
									 t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>3) | (word_cling_64_16(f01[j2]<<56)>>2) | (word_cling_64_16(f02[j2]<<56)>>1) | (word_cling_64_16(f03[j2]<<56)>>0);
									 t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>3) | (word_cling_64_16(f01[j2]<<52)>>2) | (word_cling_64_16(f02[j2]<<52)>>1) | (word_cling_64_16(f03[j2]<<52)>>0);
									 t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>3) | (word_cling_64_16(f01[j2]<<48)>>2) | (word_cling_64_16(f02[j2]<<48)>>1) | (word_cling_64_16(f03[j2]<<48)>>0);
									 t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>3) | (word_cling_64_16(f01[j2]<<44)>>2) | (word_cling_64_16(f02[j2]<<44)>>1) | (word_cling_64_16(f03[j2]<<44)>>0);
									 t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>3) | (word_cling_64_16(f01[j2]<<40)>>2) | (word_cling_64_16(f02[j2]<<40)>>1) | (word_cling_64_16(f03[j2]<<40)>>0);
									 t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>3) | (word_cling_64_16(f01[j2]<<36)>>2) | (word_cling_64_16(f02[j2]<<36)>>1) | (word_cling_64_16(f03[j2]<<36)>>0);
									 t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>3) | (word_cling_64_16(f01[j2]<<32)>>2) | (word_cling_64_16(f02[j2]<<32)>>1) | (word_cling_64_16(f03[j2]<<32)>>0);
									 t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>3) | (word_cling_64_16(f01[j2]<<28)>>2) | (word_cling_64_16(f02[j2]<<28)>>1) | (word_cling_64_16(f03[j2]<<28)>>0);
									 t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>3) | (word_cling_64_16(f01[j2]<<24)>>2) | (word_cling_64_16(f02[j2]<<24)>>1) | (word_cling_64_16(f03[j2]<<24)>>0);
									 t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>3) | (word_cling_64_16(f01[j2]<<20)>>2) | (word_cling_64_16(f02[j2]<<20)>>1) | (word_cling_64_16(f03[j2]<<20)>>0);
									 t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>3) | (word_cling_64_16(f01[j2]<<16)>>2) | (word_cling_64_16(f02[j2]<<16)>>1) | (word_cling_64_16(f03[j2]<<16)>>0);
									 t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>3) | (word_cling_64_16(f01[j2]<<12)>>2) | (word_cling_64_16(f02[j2]<<12)>>1) | (word_cling_64_16(f03[j2]<<12)>>0);
									 t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>3) | (word_cling_64_16(f01[j2]<< 8)>>2) | (word_cling_64_16(f02[j2]<< 8)>>1) | (word_cling_64_16(f03[j2]<< 8)>>0);
									 t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>3) | (word_cling_64_16(f01[j2]<< 4)>>2) | (word_cling_64_16(f02[j2]<< 4)>>1) | (word_cling_64_16(f03[j2]<< 4)>>0);
									 t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>3) | (word_cling_64_16(f01[j2]<< 0)>>2) | (word_cling_64_16(f02[j2]<< 0)>>1) | (word_cling_64_16(f03[j2]<< 0)>>0);
								 }

								 register word tmp = t[T->x->width-1];
								 switch(T->x->width - j) {
									 case 16: t[j+15] |= (word_cling_64_16(f00[j2]<< 0)>>3) | (word_cling_64_16(f01[j2]<< 0)>>2) | (word_cling_64_16(f02[j2]<< 0)>>1) | (word_cling_64_16(f03[j2]<< 0)>>0);
									 case 15: t[j+14] |= (word_cling_64_16(f00[j2]<< 4)>>3) | (word_cling_64_16(f01[j2]<< 4)>>2) | (word_cling_64_16(f02[j2]<< 4)>>1) | (word_cling_64_16(f03[j2]<< 4)>>0);
									 case 14: t[j+13] |= (word_cling_64_16(f00[j2]<< 8)>>3) | (word_cling_64_16(f01[j2]<< 8)>>2) | (word_cling_64_16(f02[j2]<< 8)>>1) | (word_cling_64_16(f03[j2]<< 8)>>0);
									 case 13: t[j+12] |= (word_cling_64_16(f00[j2]<<12)>>3) | (word_cling_64_16(f01[j2]<<12)>>2) | (word_cling_64_16(f02[j2]<<12)>>1) | (word_cling_64_16(f03[j2]<<12)>>0);
									 case 12: t[j+11] |= (word_cling_64_16(f00[j2]<<16)>>3) | (word_cling_64_16(f01[j2]<<16)>>2) | (word_cling_64_16(f02[j2]<<16)>>1) | (word_cling_64_16(f03[j2]<<16)>>0);
									 case 11: t[j+10] |= (word_cling_64_16(f00[j2]<<20)>>3) | (word_cling_64_16(f01[j2]<<20)>>2) | (word_cling_64_16(f02[j2]<<20)>>1) | (word_cling_64_16(f03[j2]<<20)>>0);
									 case 10: t[j+ 9] |= (word_cling_64_16(f00[j2]<<24)>>3) | (word_cling_64_16(f01[j2]<<24)>>2) | (word_cling_64_16(f02[j2]<<24)>>1) | (word_cling_64_16(f03[j2]<<24)>>0);
									 case  9: t[j+ 8] |= (word_cling_64_16(f00[j2]<<28)>>3) | (word_cling_64_16(f01[j2]<<28)>>2) | (word_cling_64_16(f02[j2]<<28)>>1) | (word_cling_64_16(f03[j2]<<28)>>0);
									 case  8: t[j+ 7] |= (word_cling_64_16(f00[j2]<<32)>>3) | (word_cling_64_16(f01[j2]<<32)>>2) | (word_cling_64_16(f02[j2]<<32)>>1) | (word_cling_64_16(f03[j2]<<32)>>0);
									 case  7: t[j+ 6] |= (word_cling_64_16(f00[j2]<<36)>>3) | (word_cling_64_16(f01[j2]<<36)>>2) | (word_cling_64_16(f02[j2]<<36)>>1) | (word_cling_64_16(f03[j2]<<36)>>0);
									 case  6: t[j+ 5] |= (word_cling_64_16(f00[j2]<<40)>>3) | (word_cling_64_16(f01[j2]<<40)>>2) | (word_cling_64_16(f02[j2]<<40)>>1) | (word_cling_64_16(f03[j2]<<40)>>0);
									 case  5: t[j+ 4] |= (word_cling_64_16(f00[j2]<<44)>>3) | (word_cling_64_16(f01[j2]<<44)>>2) | (word_cling_64_16(f02[j2]<<44)>>1) | (word_cling_64_16(f03[j2]<<44)>>0);
									 case  4: t[j+ 3] |= (word_cling_64_16(f00[j2]<<48)>>3) | (word_cling_64_16(f01[j2]<<48)>>2) | (word_cling_64_16(f02[j2]<<48)>>1) | (word_cling_64_16(f03[j2]<<48)>>0);
									 case  3: t[j+ 2] |= (word_cling_64_16(f00[j2]<<52)>>3) | (word_cling_64_16(f01[j2]<<52)>>2) | (word_cling_64_16(f02[j2]<<52)>>1) | (word_cling_64_16(f03[j2]<<52)>>0);
									 case  2: t[j+ 1] |= (word_cling_64_16(f00[j2]<<56)>>3) | (word_cling_64_16(f01[j2]<<56)>>2) | (word_cling_64_16(f02[j2]<<56)>>1) | (word_cling_64_16(f03[j2]<<56)>>0);
									 case  1: t[j+ 0] |= (word_cling_64_16(f00[j2]<<60)>>3) | (word_cling_64_16(f01[j2]<<60)>>2) | (word_cling_64_16(f02[j2]<<60)>>1) | (word_cling_64_16(f03[j2]<<60)>>0);
											  break;
									 default:
											  m4ri_die("impossible");
								 }
								 t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
							 }
						 }
						 break;
			}
		}
	}
	return T;
}
