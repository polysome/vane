#ifndef _GF2EXT_HPP_
#define _GF2EXT_HPP_

#include <stdint.h>
#include <stdio.h>

#include <openssl/rand.h>




template <unsigned W> struct gf2ext_u64;

template <unsigned W> inline
gf2ext_u64<W> _mul( const gf2ext_u64<W> & a, const gf2ext_u64<W> & b );





template <unsigned W>
struct gf2ext_u64 {
	uint64_t v[2];

	typedef gf2ext_u64 gfext;
	typedef gf2ext_u64 gf_t;

	static const unsigned _num_byte = (W+7)/8;

	static const uint64_t _zero[2];
	static const uint64_t _one[2];
	static const uint64_t _msb_u64 __attribute__((aligned(16))) = 0x8000000000000000ULL;

	static const uint64_t _mask_high_limb __attribute__((aligned(16))) = 0xffffffffffffffffULL>>(128-W);
	static const unsigned _nbit_high_limb = (W-64);
	static const uint64_t _irrPoly[2];

	gf2ext_u64() { v[0]^=v[0]; v[1]^=v[1]; }

	gf2ext_u64(const gf2ext_u64 & a ) { v[0]=a.v[0]; v[1]=a.v[1]; }

	gf2ext_u64( const uint8_t * x ) {
		v[0]=((const uint64_t *)x)[0];
		for(unsigned i=8;i<_num_byte;i++) ((uint8_t*)v)[i]=x[i];
		v[1] &= _mask_high_limb;
	}

	gf_t inv() const; /// XXX:
	gf_t squ() const { return _mul( *this , *this ); } /// XXX:

	gf_t & operator *= ( const gf_t & a ) { *this = _mul( *this, a ); return *this; }
	gf_t & operator ^= ( const gf_t & a ) { v[0]^=a.v[0]; v[1]^=a.v[1]; return *this; }

	gf_t operator *( const gf_t & a ) const { return _mul(*this,a); }
	gf_t operator ^( const gf_t & a ) const { gf_t r=*this; r^=a; return r; }

	gf_t & set_zero() { v[0]^=v[0]; v[1]^=v[1]; return *this; }
	bool is_zero() const { return !(v[0]|v[1]); }
	bool is_one() const { return (0==v[1])&&(1==v[0]); }

	static gf_t rand() {
		gf_t r;
		RAND_bytes( (unsigned char*)r.v, 16 );
		r.v[1] &= _mask_high_limb;
		return r;
	}
	static const gf_t& zero() { return (const gf_t&)_zero; }
	static const gf_t& one() { return (const gf_t&)_one; }
	static const gf_t& irrPoly() { return (const gf_t&)_irrPoly; }

	static gf_t assign( const uint8_t * x ) { return gf_t(x); }

	static unsigned num_byte() { return _num_byte; }
	void dump( uint8_t * x ) const { for(unsigned i=0;i<_num_byte;i++) x[i]=((uint8_t*)v)[i]; }

	void fdump(FILE *fp) const {
		uint16_t *v16 = (uint16_t *)v;
		if(v[1]) fprintf(fp,".%04x.",v16[4]);
		if(v[0]) fprintf(fp,".%04x",v16[0]);
		if( is_zero() ) fprintf(fp,"0");
	}

	void fdump2(FILE *fp) const {
		uint16_t *v16 = (uint16_t *)v;
		for(int i=7;i>=0;i--) fprintf(fp,"[%4x]",v16[i]);
//		if(v[1]) fprintf(fp,".%04x.",v16[4]);
//		if(v[0]) fprintf(fp,".%04x",v16[0]);
//		if( is_zero() ) fprintf(fp,"0");
	}

//// private

	gf2ext_u64 & mul_2() {
		v[1] <<= 1;
		v[1] |= (v[0]&_msb_u64)? 1:0;
		v[0] <<= 1;
		(*this) ^= ((v[1] & _irrPoly[1])? (const gf2ext_u64 &)_irrPoly : (const gf2ext_u64 &)_zero);
		return *this;
	}

	int left_most_bit() const {
		if(v[1]) return 64+63-__builtin_clzll(v[1]);
		if(v[0]) return 63-__builtin_clzll(v[0]);
		return -1;
	}

	gf2ext_u64 & bit_shl( unsigned i) {
		if( i >= 64 ) { v[1] = v[0]; v[0] = 0; i -= 64; }
		v[1] = (v[1]<<i)| (v[0]>>(64-i));
		v[0] <<= i;
		return *this;
	}
};



template <unsigned W>
const uint64_t gf2ext_u64<W>::_zero[2] __attribute__((aligned(16))) = {0};

template <unsigned W>
const uint64_t gf2ext_u64<W>::_one[2] __attribute__((aligned(16))) = {1ULL,0ULL};


#include "run_config.h"

#ifdef CONFIG_PROFILE
#include "stat_profile.h"
#endif


template <unsigned W> inline
gf2ext_u64<W> _mul( const gf2ext_u64<W> & a, const gf2ext_u64<W> & b )
{
#ifdef CONFIG_PROFILE
count_mul();
#endif
	gf2ext_u64<W> bpower = b;
	gf2ext_u64<W> accu;

	uint64_t a0 = a.v[0];
	for( int i=0;i<64;i++) {
		if( a0&0x01 ) accu ^= bpower;
		a0 >>= 1;
		bpower.mul_2();
	}
	a0 = a.v[1];
	for( unsigned i=0;i< gf2ext_u64<W>::_nbit_high_limb;i++) {
		if( a0&0x01 ) accu ^= bpower;
		a0 >>= 1;
		bpower.mul_2();
	}

	return accu;
}


template <unsigned W>
gf2ext_u64<W> gf2ext_u64<W>::inv() const
{
#ifdef CONFIG_PROFILE
count_inv();
#endif
	gfext buf1[2];
	gfext buf2[2];
	gfext * ptr1 = buf1;
	gfext * ptr2 = buf2;
	gfext * tmp;

	buf1[0] = one();
	buf1[1] = *this;
	int lmb1 = left_most_bit();
	buf2[0] = zero();
	buf2[1] = irrPoly();
	int lmb2 = W;
	gfext a;
	gfext b;
	while( lmb1 > 0 ) {
		a = ptr1[0];
		b = ptr1[1];
		unsigned diff = lmb2-lmb1;
		if(diff) a.bit_shl(diff);
		if(diff) b.bit_shl(diff);
		ptr2[0] ^= a;
		ptr2[1] ^= b;
		lmb2 = ptr2[1].left_most_bit();
		if( lmb2 < lmb1 ) {
			int t = lmb1; lmb1=lmb2; lmb2=t;
			tmp = ptr1; ptr1=ptr2; ptr2=tmp;
		}
	}
	return ptr1[0];
}



#include "run_config.h"


#ifdef CONFIG_HAS_PCLMULQDQ
#include "gf2ext-sse.hpp"
#define GF2EXT gf2ext_sse
#else
#define GF2EXT gf2ext_u64
#endif



#endif /// _GF2EXT_HPP_
