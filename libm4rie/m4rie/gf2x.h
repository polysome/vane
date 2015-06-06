/**
 * \file gf2x.h
 *
 * \brief \GF2X for degrees < 64
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 *
 * \warning Do not rely on these functions for high performance, they are not fully optimised.
 */

#ifndef M4RIE_GF2X_H
#define M4RIE_GF2X_H

/******************************************************************************
 *
 *            M4RIE: Linear Algebra over GF(2^e)
 *
 *    Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include <m4ri/m4ri.h>

#define __M4RIE_1tF(X) ~((X)-1) /**< Maps 1 to word with all ones and 0 to 0. */

typedef int deg_t; /**< degree type **/

/**
 * \brief a*b in \GF2X with deg(a) and deg(b) < d.
 */

static inline word gf2x_mul(const word a, const word b, deg_t d) {
	word res = 0;

	switch(d) {
		case 32: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(31)) >>31) & (b<<31);
		case 31: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(30)) >>30) & (b<<30);
		case 30: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(29)) >>29) & (b<<29);
		case 29: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(28)) >>28) & (b<<28);
		case 28: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(27)) >>27) & (b<<27);
		case 27: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(26)) >>26) & (b<<26);
		case 26: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(25)) >>25) & (b<<25);
		case 25: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(24)) >>24) & (b<<24);
		case 24: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(23)) >>23) & (b<<23);
		case 23: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(22)) >>22) & (b<<22);
		case 22: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(21)) >>21) & (b<<21);
		case 21: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(20)) >>20) & (b<<20);
		case 20: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(19)) >>19) & (b<<19);
		case 19: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(18)) >>18) & (b<<18);
		case 18: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(17)) >>17) & (b<<17);
		case 17: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(16)) >>16) & (b<<16);
		case 16: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(15)) >>15) & (b<<15);
		case 15: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(14)) >>14) & (b<<14);
		case 14: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(13)) >>13) & (b<<13);
		case 13: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(12)) >>12) & (b<<12);
		case 12: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(11)) >>11) & (b<<11);
		case 11: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW(10)) >>10) & (b<<10);
		case 10: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 9)) >> 9) & (b<< 9);
		case  9: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 8)) >> 8) & (b<< 8);
		case  8: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 7)) >> 7) & (b<< 7);
		case  7: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 6)) >> 6) & (b<< 6);
		case  6: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 5)) >> 5) & (b<< 5);
		case  5: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 4)) >> 4) & (b<< 4);
		case  4: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 3)) >> 3) & (b<< 3);
		case  3: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 2)) >> 2) & (b<< 2);
		case  2: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 1)) >> 1) & (b<< 1);
		case  1: res ^= __M4RIE_1tF((a & __M4RI_TWOPOW( 0)) >> 0) & (b<< 0);
				 break;
		default:
				 m4ri_die("degree %d too big.\n",d);
	}
	return res;
}
/**
 * \brief deg(a) in \GF2X.
 *
 * \param a Polynomial of degree <= 64.
 */

static inline deg_t gf2x_deg(word a) {
	deg_t degree = 0;
	if( (a & 0xffffffff00000000ULL) != 0) { degree += 32; a>>=32; }
	if( (a &         0xffff0000ULL) != 0) { degree += 16; a>>=16; }
	if( (a &             0xff00ULL) != 0) { degree +=  8; a>>= 8; }
	if( (a &               0xf0ULL) != 0) { degree +=  4; a>>= 4; }
	if( (a &                0xcULL) != 0) { degree +=  2; a>>= 2; }
	if( (a &                0x2ULL) != 0) { degree +=  1; a>>= 1; }
	return degree;
}

/**
 * \brief a / b in \GF2X.
 */

static inline word gf2x_div(word a, word b) {
	word res = 0;
	word mask = 0;
	const deg_t deg_b = gf2x_deg(b);
	for(deg_t deg_a = gf2x_deg(a); deg_a >= deg_b; deg_a--) {
		mask = __M4RIE_1tF(a>>deg_a);
		res |= mask & __M4RI_TWOPOW(deg_a - deg_b);
		a ^=  mask & b<<(deg_a - deg_b);
	}
	return res;
}

/**
 * \brief a mod b in \GF2X.
 */

static inline word gf2x_mod(word a, word b) {
	const deg_t deg_b = gf2x_deg(b);
	for(deg_t deg_a=gf2x_deg(a); deg_a>=deg_b; deg_a--) {
		a ^= __M4RIE_1tF(a>>deg_a) & b<<(deg_a - deg_b);
	}
	return a;
}

/**
 * \brief a / b and a mod b in \GF2X.
 */

static inline word gf2x_divmod(word a, word b, word *rem) {
	word res = 0;
	word mask = 0;
	const deg_t deg_b = gf2x_deg(b);
	for(deg_t deg_a=gf2x_deg(a); deg_a>=deg_b; deg_a--) {
		mask = __M4RIE_1tF(a>>deg_a);
		res |= mask & __M4RI_TWOPOW(deg_a - deg_b);
		a ^=  mask & b<<(deg_a - deg_b);
	}
	*rem = a;
	return res;
}


/**
 * \brief a^(-1) % b with deg(a), deg(b) <= d.
 */

static inline word gf2x_invmod(word a, word b, const deg_t d) {
	word x = 0;
	word lastx = 1;
	word y = 1;
	word lasty = 0;

	word rem;
	word quo;
	word tmp = 0;

	while (b != 0) {
		quo = gf2x_divmod(a,b, &rem);
		a = b; b = rem;
		tmp = x; x = lastx ^ gf2x_mul(quo, x, d); lastx = tmp;
		tmp = y; y = lasty ^ gf2x_mul(quo, y, d); lasty = tmp;
	}
	return lastx;
}

#endif //M4RIE_GF2X_H
