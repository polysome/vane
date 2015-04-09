#ifndef _BLAS_SSE_H_
#define _BLAS_SSE_H_


#include "blas.h"

#include "emmintrin.h"

/// assume 64 < W <= 128
template <unsigned W>
struct vec_sse {
	typedef vec_sse vec_t;
	typedef vec_u64x2<W> v_u64;
	__m128i v;

	vec_sse() { v = _mm_setzero_si128(); }
	vec_sse(const uint8_t *inp) {
		v_u64 & v64 = (v_u64 &)(*this);
		v64 = v_u64(inp);
	}
	vec_sse(const vec_sse & a ) : v(a.v) {}

	static unsigned num_byte() { return v_u64::num_byte(); }
	void dump(uint8_t *m) const { ((const v_u64 &)(*this)).dump(m); }
	void fdump(FILE * fp) const { ((const v_u64 &)(*this)).fdump(fp); }

	vec_t & operator ^= ( const vec_t & b ) { v ^= b.v; return *this; }

	bool get_ele(unsigned i) const {
		i&=127; int i16 = i&7; uint32_t mask = 1<<(i>>3);
		return mask&_mm_movemask_epi8( _mm_slli_epi16(v,7-i16) );
	}
	void set_ele(unsigned i) { ((v_u64 &)(*this)).set_ele(i); }

	template <unsigned tail_bit>
	uint64_t tail () const { return ((const v_u64 &)(*this)). template tail<tail_bit>(); }
	template <unsigned W2>
	vec_sse<W2> concate( uint64_t tail = 0 ) const { vec_sse<W2> r; *(vec_u64x2<W2>*)(&r) = ((const v_u64&)(*this)). template concate<W2>( tail ); return r; }

	vec_t& set_zero() { v = _mm_setzero_si128(); return *this; }
	bool is_zero() const { return (0xffff==_mm_movemask_epi8(_mm_cmpeq_epi16(_mm_setzero_si128(),v))); }
	bool operator==( const vec_t & a ) const { return (0xffff==_mm_movemask_epi8(_mm_cmpeq_epi16(a.v,v))); }

	static vec_t rand() { vec_t r; *(v_u64*)(&r) = v_u64::rand(); return r; }
	static const vec_t& zero() { return (const vec_t &)vec_u64x4_zero; }

};


#endif
