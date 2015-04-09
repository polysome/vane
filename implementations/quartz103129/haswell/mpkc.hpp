
#ifndef _MPKC_HPP_
#define _MPKC_HPP_

#include "blas.h"



#define TERMS_QUAD_POLY(N) ((N)*(N+1)/2)

#define IDX_XSQ(i,n_var) ((2*(n_var)+1-i)*(i)/2)


template <unsigned n_var,unsigned n_eqs>
void interpolate( MAT<TERMS_QUAD_POLY(n_var),n_eqs> & poly , void (*quad_poly)(void *,const void *,const void *) , const void * key )
{
	VEC<n_var> tmp;
	VEC<n_eqs> tmp_r0;
	VEC<n_eqs> tmp_r1;

	for(unsigned i=0;i<n_var;i++) {
		unsigned base_idx = IDX_XSQ(i,n_var);
		tmp.set_zero();
		tmp.set_ele(i);
		quad_poly( &tmp_r0 , key , & tmp );
		poly.v[base_idx] = tmp_r0;
	}

	for(unsigned i=0;i<n_var;i++) {
		unsigned base_idx = IDX_XSQ(i,n_var);
		tmp_r0 = poly.v[base_idx];

		for(unsigned j=i+1;j<n_var;j++) {
			tmp.set_zero();
			tmp.set_ele(i);
			tmp.set_ele(j);

			quad_poly( &tmp_r1 , key , & tmp );
			tmp_r1 ^= tmp_r0;
			tmp_r1 ^= poly.v[IDX_XSQ(j,n_var)];

			poly.v[base_idx+j-i] = tmp_r1;
		}
	}
}


template <unsigned n_var,unsigned n_eqs>
void mpkc_pub_map( VEC<n_eqs> & z , const MAT<TERMS_QUAD_POLY(n_var),n_eqs> & pk , const VEC<n_var> & w )
{
	unsigned acc_idx = 0;
	VEC<n_eqs> r;

	for(unsigned i=0;i<n_var;i++) {
		if( ! w.get_ele(i) ) { acc_idx += (n_var-i); continue; }
		r ^= pk.v[acc_idx++];
		for(unsigned j=i+1;j<n_var;j++,acc_idx++) if( w.get_ele(j) ) r^= pk.v[acc_idx];
	}
	z = r;
}



#endif
