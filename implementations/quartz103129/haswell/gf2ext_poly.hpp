#ifndef _GF2EXT_POLY_HPP_
#define _GF2EXT_POLY_HPP_


#include "gf2ext.hpp"



#include "gf2ext.hpp"

#ifdef CONFIG_PROFILE
#include "stat_profile.h"
#endif


template <unsigned ext,unsigned max_deg>
struct poly
{
	static const unsigned max_term = max_deg+1;

	typedef GF2EXT<ext> gf_t;
	typedef poly<ext,max_deg> poly_t;

	gf_t coef[max_term];
	unsigned deg[max_term];
	unsigned n_terms;

	poly() : n_terms(0) {}

	void set_zero() { n_terms = 0; }

	bool append_term( const gf_t & c , unsigned d ) {
		if(max_term == n_terms ) return false;
		coef[n_terms]=c;
		deg[n_terms]=d;
		n_terms++;
		return true;
	}


	void X_2ext_X( poly & p ) const; /// specialize in .cpp file

	bool find_unique_root( gf_t & root ) const;

	poly & normalize() {
		if( !n_terms ) return *this;
		gf_t inv = coef[0].inv(); /// XXX
		coef[0] = gf_t::one();
		for(unsigned i=1;i<n_terms;i++) coef[i] *= inv;
//			coef[i] = mul( inv , coef[i] ); /// XXX
		return *this;
	}

	gf_t eval( const gf_t & b ) const {
		gf_t r;
		if(!n_terms) return r;

		int i = n_terms-1;
		if(0==deg[i]) { r=coef[i]; i--; }

		gf_t accu = b;
		unsigned accu_deg = 1;
		for(;i>=0;i--) {
			while( deg[i] > accu_deg ) { accu *= b; accu_deg++; }
			gf_t tt = coef[i];
			tt *= accu;
			r ^= tt;
		}
		return r;
	}

	poly & operator += ( const gf_t & c ) {
		if( (!n_terms) || (0 != deg[n_terms-1]) ) append_term( c , 0 );
		else coef[n_terms-1] ^= c;
		return *this;
	}

	poly & operator *= ( const gf_t & m ) {
		for(unsigned i=0;i<n_terms;i++) coef[i] *= m; //coef[i] = mul( m , coef[i] ); /// XXX
		return *this;
	}

	poly & mul_X(unsigned d=1) {
		for(unsigned i=0;i<n_terms;i++) deg[i]+=d;
		return *this;
	}

	poly & squ() {
		for(unsigned i=0;i<n_terms;i++) {
			deg[i] += deg[i];
			coef[i] = coef[i].squ(); /// XXX
		}
		return *this;
	}

	poly & pow2to(unsigned q) {
		for(unsigned i=0;i<n_terms;i++) {
			deg[i]<<=q;
			for(unsigned j=0;j<q;j++)
				coef[i] = coef[i].squ(); /// XXX
		}
		return *this;
	}

	void split( poly & p1 , poly & p2 , unsigned p2_max_deg ) const {
		p1.n_terms = 0;
		p2.n_terms = 0;
		for(unsigned i=0;i<n_terms;i++)
			if(deg[i]>p2_max_deg) p1.append_term(coef[i],deg[i]);
			else p2.append_term(coef[i],deg[i]);
	}

	poly & operator= ( const poly & poly ) {
		n_terms = poly.n_terms;
		for(unsigned i=0;i<n_terms;i++) {
			coef[i] = poly.coef[i];
			deg[i] = poly.deg[i];
		}
		return *this;
	}

	poly & prune() {
		unsigned new_terms = 0;
		for(unsigned i=0;i<n_terms;i++) {
			if( coef[i].is_zero() ) continue;
			if( new_terms == i ) { new_terms++; continue; }
			coef[new_terms] = coef[i];
			deg[new_terms] = deg[i];
			new_terms++;
		}
		n_terms = new_terms;
		return (*this);
	}

	poly & prune_head() {
		if( 0 == n_terms ) return *this;
		if(  coef[0].is_zero() ) return prune();
		return *this;
	}

/////////
	poly & reduce_term_by( const poly & nor_poly , poly & res )
	{
		if(!nor_poly.n_terms) return (res = *this);
		unsigned idx_a = 0;
		for (; idx_a < n_terms ;idx_a++) if( deg[idx_a] == nor_poly.deg[0] ) break;
		if( idx_a == n_terms ) return (res = *this);
		_add<true,false,true>( res , *this , nor_poly , 0 , coef[idx_a] );
		return res;
	}

	poly & reduce_headterm_by( const poly & nor_poly , poly & res )
	{
		if(!nor_poly.n_terms) return (res = *this);
		if( !n_terms ) return (res=*this);
		if( nor_poly.deg[0] > deg[0] ) return (res=*this);

		_add<true,true,true>( res , *this ,nor_poly, deg[0]-nor_poly.deg[0] , coef[0] );
		return res;
	}

	poly & reduce_headterm_add( poly & res , const poly & poly_b )
	{
		if(!poly_b.n_terms) return (res = *this);
		if( !n_terms ) return (res=*this);
		if( poly_b.deg[0] > deg[0] ) return (res=*this);

		_add<true,true,false>(res,*this,poly_b,deg[0]-poly_b.deg[0],gf_t::one());
		return res;
	}

	static void add( poly & c , const poly & a , const poly & b , const gf_t & mb )
	{
		return _add<false,false,true>(c,a,b,0,mb);
	}
////////////

	static void _euclid_gcd( poly & gcd , const poly & p1 , const poly & p2 );

	static void euclid_gcd( poly & gcd , const poly & p1 , const poly & p2 ) {
#ifdef CONFIG_PROFILE
new_phase();
BENCHMARK( ext_gcd_poly ,{
#endif
		_euclid_gcd( gcd , p1 , p2 );
#ifdef CONFIG_PROFILE
});
new_phase();
#endif
	}

	void fdump(FILE *fp) const {
		fprintf(fp,"[[%d]]",n_terms);
		if(!n_terms) return;
		fprintf(fp,"[");
		coef[0].fdump(fp);
		fprintf(fp," X^%d]",deg[0]);
		if(1==n_terms) return;
		if(n_terms>2) fprintf(fp," + ... ");
		fprintf(fp,"+ [");
		coef[n_terms-1].fdump(fp);
		if(deg[n_terms-1]) fprintf(fp," X^%d]",deg[n_terms-1]);
		else fprintf(fp,"]");
	}

	void fdump2(FILE *fp) const {
		fprintf(fp,"[[%d]]",n_terms);
		if(!n_terms) return;
		fprintf(fp,"[");
		coef[0].fdump(fp);
		fprintf(fp," X^%d]",deg[0]);
		for( unsigned i=1;i<n_terms;i++) {
			fprintf(fp," + [");
			coef[i].fdump(fp);
			if(deg[i]) fprintf(fp," X^%d]",deg[i]);
			else fprintf(fp,"]");
		}	
	}

//////////////////////////////

	template <bool annaliate_b_head, bool deg_mod, bool coef_mod >
	static void _add( poly & c , const poly & a , const poly & b , unsigned db , const gf_t & mb );

};



///////////////////////////


template <unsigned ext,unsigned max_deg,bool flag>
struct mod {
	typedef typename poly<ext,max_deg>::gf_t gf_t;

	static unsigned _deg( unsigned deg_b , unsigned deg_diff_b ) { return deg_b + deg_diff_b; }
	static gf_t _coef( const gf_t &b , const gf_t & mb ) { return b*mb; }
};
template <unsigned ext,unsigned max_deg>
struct mod<ext,max_deg,false> {
	typedef typename poly<ext,max_deg>::gf_t gf_t;

	static unsigned _deg( unsigned deg_b , unsigned ) { return deg_b; }
	static gf_t _coef( const gf_t &b , const gf_t & ) { return b; }
};


//#define _LESS_CHECK_ZERO_

template <unsigned ext,unsigned max_deg>
template <bool annaliate_b_head, bool deg_mod, bool coef_mod >
void poly<ext,max_deg>::_add( poly & c , const poly & a , const poly & b , unsigned db , const gf_t & mb ) /// XXX
{
#ifdef CONFIG_PROFILE
	count_poly_add();
#endif
	unsigned idx_a = 0;
	unsigned idx_b = 0;
	c.n_terms = 0;
	for(;idx_a<a.n_terms && idx_b<b.n_terms;) {
		if(a.deg[idx_a]== mod<ext,max_deg,deg_mod>::_deg(b.deg[idx_b],db) ) {
			if(annaliate_b_head && (0==idx_b)) { idx_a++; idx_b++; continue; }
			c.coef[c.n_terms] = a.coef[idx_a]^
					mod<ext,max_deg,coef_mod>::_coef(b.coef[idx_b],mb);
			c.deg[c.n_terms] = a.deg[idx_a];
			idx_a++;
			idx_b++;
#ifdef _LESS_CHECK_ZERO_
			c.n_terms++;
#else
			if(!c.coef[c.n_terms].is_zero()) c.n_terms++;
#endif
		} else if( a.deg[idx_a] > mod<ext,max_deg,deg_mod>::_deg(b.deg[idx_b],db) ) {
			c.coef[c.n_terms] = a.coef[idx_a];
			c.deg[c.n_terms] = a.deg[idx_a];
			c.n_terms++;
			idx_a++;
		} else {
			c.coef[c.n_terms] = mod<ext,max_deg,coef_mod>::_coef(b.coef[idx_b],mb);
			c.deg[c.n_terms] = mod<ext,max_deg,deg_mod>::_deg(b.deg[idx_b],db);
			c.n_terms++;
			idx_b++;
		}
	}
	for( ; idx_a < a.n_terms ; idx_a++) {
		c.coef[c.n_terms] = a.coef[idx_a];
		c.deg[c.n_terms] = a.deg[idx_a];
		c.n_terms++;
	}
	for( ; idx_b < b.n_terms ; idx_b++ ) {
		c.coef[c.n_terms] = mod<ext,max_deg,coef_mod>::_coef(b.coef[idx_b],mb);
		c.deg[c.n_terms] = mod<ext,max_deg,deg_mod>::_deg(b.deg[idx_b],db);
		c.n_terms++;
	}

}


//////////////////////////////////////////

#include <algorithm>


template <unsigned ext,unsigned max_deg>
void poly<ext,max_deg>::_euclid_gcd( poly<ext,max_deg> & gcd , const poly<ext,max_deg> & p1 , const poly<ext,max_deg> & p2 )
{
	if( ! p1.n_terms ) { gcd = p2; return; }
	if( ! p2.n_terms ) { gcd = p1; return; }
	poly buf1;
	poly buf2;
	poly buf3;
	buf1 = p1;
	buf2 = p2;
	poly * ptr1 = &buf1;
	poly * ptr2 = &buf2;
	poly * ptr_res = & buf3;
	poly * tmp;

	gf_t m1;
	gf_t m2;
#ifdef _LESS_CHECK_ZERO_
	ptr1->prune_head();
	ptr2->prune_head();
#endif
	for(;;){
		if(ptr1->deg[0] < ptr2->deg[0]) {
			tmp = ptr2; ptr2=ptr1; ptr1=tmp;
		}
		m1 = ptr1->coef[0];
		m2 = ptr2->coef[0];
		(*ptr1) *= m2;
		(*ptr2) *= m1;  /// XXX
		ptr1->reduce_headterm_add(*ptr_res,*ptr2);
#ifdef _LESS_CHECK_ZERO_
		ptr_res->prune();
#endif
		if( 0 == ptr_res->n_terms ) {
			gcd = *ptr2;
			return;
		}
		tmp=ptr1; ptr1=ptr_res; ptr_res=tmp;
	}
}




template <unsigned ext,unsigned max_deg>
bool poly<ext,max_deg>::find_unique_root( gf_t & root ) const
{
	poly x2103_x;

	X_2ext_X( x2103_x );

	poly gcd;

#ifdef CONFIG_PROFILE
	BENCHMARK( ext_gcd_poly , {
#endif

	euclid_gcd( gcd , *this , x2103_x );

#ifdef CONFIG_PROFILE
	} );
	profile_end();
#endif

	bool uni_root = true;

	if( !gcd.n_terms ) uni_root =false;
	if( 1 != gcd.deg[0] ) uni_root = false;

	gf_t gf_one = gcd.coef[0].inv();
	root = gcd.coef[1] * gf_one;
	if( 1 == gcd.n_terms ) root = gf_t::zero();

	return uni_root;
}





#endif   /// _GFEXT_POLY_HPP_

