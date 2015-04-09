#ifndef _BLAS_H_
#define _BLAS_H_

#include <stdint.h>
#include <stdio.h>
#include <openssl/rand.h>

#include "run_config.h"

extern uint64_t vec_u64x4_zero[4];

/// assume 64 < W <= 128
template <unsigned W>
struct vec_u64x2 {
	typedef vec_u64x2 vec_t;
	static const unsigned _num_byte = (W+7)/8;
	static const unsigned _num_limb = (W+63)/64;
	static const unsigned _idx_last_limb = (W-1)/64;
	static const uint64_t _high_mask = (0==(W&63))? (~0ULL) : ((0x1ULL)<<(W&63))-1;
	uint64_t v[_num_limb];

	vec_u64x2() {for(unsigned i=0;i<_num_limb;i++)v[i]=0;}
	vec_u64x2(const uint8_t *inp) {
		for(unsigned i=0;i<_num_limb-1;i++)v[i]=((const uint64_t *)inp)[i];
		for(unsigned i=_idx_last_limb*8;i<_num_byte;i++) ((uint8_t*)v)[i]=inp[i];
		v[_idx_last_limb] &= _high_mask;
	}
	vec_u64x2( const vec_t & a ) {for(unsigned i=0;i<_num_limb;i++)v[i]=a.v[i];}

	static unsigned num_byte() { return _num_byte; }
	void dump(uint8_t *m) const { for(unsigned i=0;i<_num_byte;i++) m[i]=((uint8_t*)v)[i]; }
	void fdump(FILE * fp) const { uint32_t *v32=(uint32_t*)v;
		fprintf(fp,"[%3d][",W);
		for(int i=_idx_last_limb;i>=0;i--) fprintf(fp,"0x%x,0x%x,",v32[i*2+1],v32[i*2]);
		fprintf(fp,"]");
	}

	vec_t & operator ^= ( const vec_t & b ) {for(unsigned i=0;i<_num_limb;i++)v[i]^=b.v[i]; return *this;}

	bool get_ele(unsigned i) const {
		if(i>=W) return false;
		uint64_t mask = 1ULL<<(i&63);
		return (v[i>>6]&mask)==mask;
	}
	void set_ele(unsigned i) {
		if(i>=W) return;
		uint64_t mask = 1ULL<<(i&63);
		v[i>>6] |= mask;
	}

	template <unsigned tail_bit>
	uint64_t tail () const {
		if(W>128) { uint64_t t=(v[2]<<16)|(v[1]>>48); return t>>(W-112-tail_bit); }
		if(W>64) return v[1]>>(W-64-tail_bit);
		printf("unsupported tail()\n");
		exit(-1);
		return 0;
	}

	template <unsigned W2>
	vec_u64x2<W2> concate( uint64_t tail = 0 ) const {
		if( W<=128 && W2 > 128 ) {
			vec_u64x2<W2> r;
			r.v[0]=v[0];
			r.v[1]=v[1] | tail<<(64-(128-W));
			r.v[2]=tail>>(128-W);
			r.v[2] &= vec_u64x2<W2>::_high_mask;
			return r;
		}
		typedef vec_u64x2<W2> tar_t;
		tar_t r((const uint8_t *)this);
		uint64_t t= tail<<(W&63);
		r.v[tar_t::_idx_last_limb] |= t;
		r.v[tar_t::_idx_last_limb] &= tar_t::_high_mask;
		return r;
	}

	vec_t& set_zero() { *this^=*this; return *this; }

	bool is_zero() const { if(W>128) return !(v[0]|v[1]|v[2]); return !(v[0]|v[1]); }
	bool operator==( const vec_t & a ) const { vec_t t=a; t^=*this; return t.is_zero(); }

	static vec_t rand() {
		vec_t r;
		RAND_bytes( (unsigned char*)r.v , _num_byte );
		r.v[_idx_last_limb] &= _high_mask;
		return r;
	}
	static const vec_t& zero() { return (const vec_t &)vec_u64x4_zero; }
};





#ifdef CONFIG_SSE_VEC
#include "blas-sse.h"
#define VEC vec_sse
#else
#define VEC vec_u64x2
#endif



template <unsigned N_VEC,unsigned L_VEC>
struct mat_u64{
	typedef mat_u64 mat_t;
	typedef VEC<L_VEC> vec_v;
	typedef VEC<N_VEC> vec_h;
	vec_v v[N_VEC];

	void fdump(FILE *fp) const { for(unsigned i=0;i<N_VEC;i++) {v[i].fdump(fp); fprintf(fp,"\n"); } }
	static uint32_t num_byte() { return N_VEC*vec_v::num_byte(); }
	void dump(uint8_t * m) const { for(unsigned i=0;i<N_VEC;i++) { v[i].dump(m); m += vec_v::num_byte(); } }
	void set(const uint8_t * m) { for(unsigned i=0;i<N_VEC;i++) { v[i]=vec_v(m); m += vec_v::num_byte(); } }

	mat_t & set_zero() { for(unsigned i=0;i<N_VEC;i++) v[i].set_zero(); return *this; }
	mat_t & operator ^= (const mat_t & a) { for(unsigned i=0;i<N_VEC;i++) v[i]^=a.v[i]; return *this; }

	vec_v prod( const vec_h & a_ ) const {
		vec_v r;
		r.set_zero();
		for(unsigned i=0;i<N_VEC;i++) r ^= a_.get_ele(i)?v[i]:vec_v::zero();
		return r;
	}

	static bool rand_inv( mat_u64<L_VEC,L_VEC> & a , mat_u64<L_VEC,L_VEC> & b );

	static void mul( mat_u64<L_VEC,L_VEC> & c , const mat_u64<N_VEC,L_VEC> & a , const mat_u64<L_VEC,N_VEC> & b ) {
		for(unsigned k=0;k<L_VEC;k++){
			c.v[k].set_zero();
			for(unsigned i=0;i<N_VEC;i++) c.v[k] ^= b.v[k].get_ele(i)?a.v[i]:vec_v::zero();
		}
	}
};


template <unsigned W>
bool gauss_elim( mat_u64<W,W> & a , mat_u64<W,W> & b )
{
	for(unsigned i=0;i<W;i++) {
		if(!a.v[i].get_ele(i)) {
			for(unsigned j=i+1;j<W;j++) {
				if(!a.v[j].get_ele(i)) continue;
				typename mat_u64<W,W>::vec_h tt = a.v[i]; a.v[i]=a.v[j]; a.v[j]=tt;
				tt = b.v[i]; b.v[i]=b.v[j]; b.v[j]=tt;
				break;
			}
			if(!a.v[i].get_ele(i)) return false;
		}
		for(unsigned j=0;j<W;j++) {
			if(i==j) continue;
			if(!a.v[j].get_ele(i)) continue;
			a.v[j] ^= a.v[i];
			b.v[j] ^= b.v[i];
		}
	}

	return true;
}

template <unsigned W,unsigned H>
bool mat_u64<W,H>::rand_inv( mat_u64<H,H> & a , mat_u64<H,H> & b )
{
	mat_u64<H,H> aa;
	unsigned k;
	for(k=0;k<100;k++){
		b.set_zero();
		for(unsigned i=0;i<H;i++){
			a.v[i] = mat_u64<H,H>::vec_h::rand();
			b.v[i].set_ele(i);
		}
		aa = a;
		if( gauss_elim(aa,b) ) break;
	}
	return (100!=k);
}




#define MAT mat_u64


#endif
