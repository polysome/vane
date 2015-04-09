#ifndef _GF2EXT_SSE_HPP_
#define _GF2EXT_SSE_HPP_

#include <stdint.h>


#include "emmintrin.h"
#include "wmmintrin.h"



#ifdef CONFIG_PROFILE
#include "stat_profile.h"
#endif

#include "gf2ext.hpp"



template <unsigned W> struct gf2ext_sse;

template <unsigned W> inline
gf2ext_sse<W> _mul( const gf2ext_sse<W> & a , const gf2ext_sse<W> & b );




template <unsigned W>
struct gf2ext_sse {
	__m128i v;

	typedef gf2ext_sse gf_t;
	typedef gf2ext_u64<W> gf_u64;

	gf2ext_sse() { v^=v; }
	gf2ext_sse( __m128i a ): v(a) {}
	gf2ext_sse( const gf_t & a ): v(a.v) {}
	gf2ext_sse( const gf_u64 & a ) { *((gf_u64*)this)=a; }
	gf2ext_sse( const uint8_t * x ) { *((gf_u64*)this)=gf_u64(x); }

	gf_t squ() const;
	gf_t inv() const {
		gf_t r;
		r=((const gf_u64 *)(this))->inv();
		return r;
	}

	gf_t & operator *= ( const gf_t & a ) { *this = _mul( *this, a ); return *this; }
	gf_t & operator ^= ( const gf_t & a ) { v^=a.v; return *this; }

	gf_t operator *( const gf_t & a ) const { return _mul(*this,a); }
	gf_t operator ^( const gf_t & a ) const { gf_t r=*this; r^=a; return r; }

	gf_t & set_zero() { v^=v; return *this; }
	bool is_zero() const { return ((const gf_u64 *)this)->is_zero(); }
	bool is_one() const { return ((const gf_u64 *)this)->is_one(); }

	static gf_t rand() { gf_t r = gf_u64::rand(); return r; }

	static const gf_t & zero() { return (const gf_t &)gf_u64::_zero; }
	static const gf_t & one() { return (const gf_t &)gf_u64::_one; }
	static const gf_t & irrPoly() { return (const gf_t &)gf_u64::_irrPoly; }

	static gf_t assign( const uint8_t * x ) { return gf_t(x); }
	static unsigned num_byte() { return gf_u64::num_byte(); }
	void dump( uint8_t * x ) const { ((const gf_u64 *)this)->dump(x); }
	void fdump(FILE *fp) const { ((const gf_u64 *)this)->fdump(fp); }
	void fdump2(FILE *fp) const { ((const gf_u64 *)this)->fdump2(fp); }

//	static const uint32_t _mask_56bit[4];

	static const uint32_t _mask_Wbit[4];
	static const uint64_t _reducer_W[2];


};



//__m128i _mm_clmulepi64_si128(__m128i v1, __m128i v2, const int imm8);


////////////////////////////////////////////////////////

template <unsigned W> inline 
void _reduce( __m128i &a0b0 , __m128i a1b1 );



template <unsigned W> inline
gf2ext_sse<W> gf2ext_sse<W>::squ() const
{
#ifdef CONFIG_PROFILE
	count_squ();
#endif
	__m128i a0b0 = _mm_clmulepi64_si128( v , v , 0x00 );
	__m128i a1b1 = _mm_clmulepi64_si128( v , v , 0x11 );

	_reduce<W>( a0b0 , a1b1 );

	gf_t r;
	r.v = a0b0;
	return r;
}



template <unsigned W> inline
gf2ext_sse<W> _mul( const gf2ext_sse<W> & a, const gf2ext_sse<W> & b )
{
#ifdef CONFIG_PROFILE
	count_mul();
#endif
	__m128i a0b0 = _mm_clmulepi64_si128( a.v , b.v , 0x00 );
	__m128i a1b1 = _mm_clmulepi64_si128( a.v , b.v , 0x11 );
// cross terms
#ifdef CONFIG_FAST_PCLMULQDQ
	__m128i a0b1 = _mm_clmulepi64_si128( a.v , b.v , 0x10 ) ^ _mm_clmulepi64_si128( a.v , b.v , 0x01 );
#else
	__m128i a0b1 = (__m128i)_mm_shuffle_pd( (__m128d)a.v,(__m128d)b.v,1);  /// --> [b_low64,a_high64]
	a0b1 = _mm_clmulepi64_si128( a0b1^a.v,a0b1^b.v, 0x10 ) ^ a0b0 ^ a1b1;
#endif
	a0b0 ^= (__m128i)_mm_shuffle_pd( (__m128d)_mm_setzero_si128() , (__m128d)a0b1 , 0 ); /// --> [a0b1_low,0]
	a1b1 ^= (__m128i)_mm_shuffle_pd( (__m128d)a0b1 , (__m128d)_mm_setzero_si128() , 1 ); /// --> [0,a0b0_high]
// end cross terms

	_reduce<W>( a0b0 , a1b1 );

	gf2ext_sse<W> r;
	r.v = a0b0;
	return r;
}




////////////////////////////    gf2103   reduce   ////////////////////////////////////////



inline void reduce_clmul_103( __m128i &a0b0 , __m128i a1b1 )  /// test fail. check!!!
{
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<103>::_mask_Wbit );
	__m128i rder = _mm_load_si128( (__m128i*) gf2ext_sse<103>::_reducer_W );

	a0b0 ^= _mm_clmulepi64_si128( a1b1 , rder , 0x00 )
		^ _mm_slli_si128( _mm_clmulepi64_si128( a1b1 , rder , 0x01) , 8);
	a0b0 ^= _mm_clmulepi64_si128( _mm_srli_epi64(a0b0,39) , rder , 0x11 );
	//a0b0 ^= _mm_clmulepi64_si128( _mm_andnot_si128(msk,a0b0) , rder , 0x11 );
	a0b0 &= msk;

}


inline void reduce_shift( __m128i &a0b0 , __m128i a1b1 )
{
	//__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<103>::_mask_56bit );
	// {0xffffffff,0x00ffffff,0,0};
	__m128i msk = _mm_set_epi32( 0 , 0 , 0x00ffffff , 0xffffffff );
	__m128i a1b1h = _mm_andnot_si128( msk , a1b1 );
	a1b1 &= msk;

	__m128i shr103 = _mm_slli_epi64(a1b1,1);
	__m128i shr94 = _mm_slli_epi64(a1b1,2);
	shr103 = _mm_slli_si128(shr103,3);
	shr94 = _mm_slli_si128(shr94,4);
	a0b0 ^= shr103;
	a0b0 ^= shr94;
	shr103 = _mm_slli_si128(a1b1h,3);
	shr94 = _mm_slli_si128(a1b1h,4);
	shr103 = _mm_slli_epi64(shr103,1);
	shr94 = _mm_slli_epi64(shr94,2);
	a0b0 ^= shr103;
	a0b0 ^= shr94;

	msk = _mm_load_si128( (__m128i*) gf2ext_sse<103>::_mask_Wbit );
	__m128i a0b0h = _mm_andnot_si128( msk , a0b0 );
	a0b0 &= msk;
	a0b0h = _mm_srli_si128(a0b0h,12);
	a0b0 ^= _mm_srli_epi64(a0b0h,7);
	a0b0 ^= _mm_slli_epi64(a0b0h,2);
}


template <> inline 
void _reduce<103>( __m128i &a0b0 , __m128i a1b1 )
{
#ifdef CONFIG_PCLMULQDQ_REDUCE
	reduce_clmul_103( a0b0 , a1b1 );
#else
	reduce_shift( a0b0 , a1b1 );
#endif
}






#endif  /// _GF2EXT_SSE_HPP_s
