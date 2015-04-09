

#include "gf2ext_poly.hpp"



typedef poly<103,129> poly129;



/////////////////

/// for polynomail mod operation

struct poly_mod_tab
{
	poly129 x129;
	poly129 x130_plus_2x[64];

	poly_mod_tab( const poly129 & p );

	poly129 & squ_mod( poly129 & p ) const;
};

//////////////////////////////////


struct poly_pow8_mod_tab
{
	poly129 x129;
	poly129 x136_plus_8x[112];

	poly_pow8_mod_tab( const poly129 & p );
	poly129 & pow8_mod( poly129 & p ) const;
	poly129 & Xpow2to103( poly129 & buf ) const;
};


struct poly_pow16_mod_tab
{
	poly129 x129;
	poly129 x144_plus_16x[120];

	poly_pow16_mod_tab( const poly129 & p );
	poly129 & pow16_mod( poly129 & p ) const;
	poly129 & Xpow2to103( poly129 & buf ) const;
};


struct poly_pow32_mod_tab
{
	poly129 x129;
	poly129 x130_plus_2x[64];
	poly129 x160_plus_32x[124];

	poly_pow32_mod_tab( const poly129 & p );
	poly129 & squ_mod( poly129 & p ) const;
	poly129 & pow32_mod( poly129 & p ) const;
	poly129 & Xpow2to103( poly129 & res ) const;
};



#ifdef CONFIG_PROFILE
#include "stat_profile.h"
#endif


//#define MODPOWER poly_pow8_mod_tab
#define MODPOWER poly_pow16_mod_tab
//#define MODPOWER poly_pow32_mod_tab


////////////////////////////////

#ifdef CONFIG_HAS_PCLMULQDQ
#define USE_OLD
#endif

#ifdef USE_OLD

template <>
void poly<103,129>::X_2ext_X( poly129 & x2103_x ) const
{
	poly129 x2103;

#ifdef CONFIG_PROFILE
	profile_start();
	bm_start(&prepare_poly);
#endif
	MODPOWER mod_power( *this );
#ifdef CONFIG_PROFILE
	bm_stop(&prepare_poly);
	new_phase();
	BENCHMARK( do_poly , {
#endif
	mod_power.Xpow2to103( x2103 );
#ifdef CONFIG_PROFILE
	} );
	new_phase();
#endif
	gf_t gf_one = gf_t::one();
	poly129 x;
	x.n_terms = 1;
	x.deg[0] = 1;
	x.coef[0] = gf_one;

	poly129::add( x2103_x , x2103 , x , gf_one );
}

#endif



////////////////////////////////////////////////////////////////////////


poly_mod_tab::poly_mod_tab( const poly129 & p )
{
	x129 = p;
	x129.normalize();

	poly129 x129_;
	x129_.n_terms = 1;
	x129_.deg[0] = 129;
	x129_.coef[0] = poly129::gf_t::one();

	poly129 temp;
//	poly129::add( temp , x129 , x129_ , poly129::gf_t::one() );
	x129.reduce_headterm_add( temp , x129_ );
	temp.mul_X();
	temp.reduce_term_by( x129 , x130_plus_2x[0] );
	for( unsigned i=1;i<64;i++ ) {
		x130_plus_2x[i] = x130_plus_2x[i-1];
		x130_plus_2x[i].mul_X();
		x130_plus_2x[i].reduce_term_by( x129 , temp );
		temp.mul_X();
		temp.reduce_term_by( x129 , x130_plus_2x[i] );
	}
}


poly129 & poly_mod_tab::squ_mod( poly129 & p ) const
{
	poly129 exceed_terms;
	poly129 buf1;
	poly129 buf2;
	poly129 * ptr1 = &buf1;
	poly129 * ptr2 = &buf2;

	p.squ();
	p.split( exceed_terms , buf1 , 128 );

	for( unsigned i=0; i<exceed_terms.n_terms; i++ ) {
		poly129::add( *ptr2 , *ptr1 , x130_plus_2x[(exceed_terms.deg[i]-130)>>1] , exceed_terms.coef[i] );
		poly129 * tmp = ptr2;
		ptr2 = ptr1;
		ptr1 = tmp;
	}
	p = *ptr1;
	return p;
}



////////////////////////////////


poly_pow8_mod_tab::poly_pow8_mod_tab( const poly129 & p )
{
	x129 = p;
	x129.normalize();

	poly129 x129_;
	x129_.n_terms = 1;
	x129_.deg[0] = 129;
	x129_.coef[0] = poly129::gf_t::one();

	poly129 buf1;
	poly129 buf2;
	poly129 * ptr1 = & buf1;
	poly129 * ptr2 = & buf2;
	poly129 * ptr_tmp;
	///x136_plus_8x[112];
	//*ptr1 = x129;
	//poly129::add( *ptr1 , x129 , x129_ , poly129::gf_t::one() );
	x129.reduce_headterm_add( *ptr1 , x129_ );
	for(int i=7;i>0;i--) {
		ptr1->mul_X();
		ptr1->reduce_term_by( x129 , *ptr2 );
		ptr_tmp=ptr1; ptr1=ptr2; ptr2=ptr_tmp;
	}
	x136_plus_8x[0] = *ptr1;
	for( unsigned i=1;i<112;i++ ) {
		for(int j=8;j>0;j--) {
			ptr1->mul_X();
			ptr1->reduce_term_by( x129 , *ptr2 );
			ptr_tmp=ptr1; ptr1=ptr2; ptr2=ptr_tmp;
		}
		x136_plus_8x[i] = *ptr1;
	}
}

poly129 & poly_pow8_mod_tab::pow8_mod( poly129 & p ) const
{
	poly129 exceed_terms;
	poly129 buf1;
	poly129 buf2;
	poly129 * ptr1 = &buf1;
	poly129 * ptr2 = &buf2;

	p.pow2to(3);
	p.split( exceed_terms , buf1 , 128 );

	for( unsigned i=0; i<exceed_terms.n_terms; i++ ) {
		poly129::add( *ptr2 , *ptr1 , x136_plus_8x[(exceed_terms.deg[i]-136)>>3] , exceed_terms.coef[i] );
		poly129 * tmp = ptr2;
		ptr2 = ptr1;
		ptr1 = tmp;
	}
	p = *ptr1;
	return p;
}


poly129 & poly_pow8_mod_tab::Xpow2to103( poly129 & x2103 ) const
{
	x2103 = x136_plus_8x[111];  // x^2^(7+3)
	for(unsigned i=10;i<103;i+=3) pow8_mod( x2103 );
	return x2103;
}


///////////////////////////////////////////


poly_pow16_mod_tab::poly_pow16_mod_tab( const poly129 & p )
{
	x129 = p;
	x129.normalize();

	poly129 x129_;
	x129_.n_terms = 1;
	x129_.deg[0] = 129;
	x129_.coef[0] = poly129::gf_t::one();

	poly129 buf1;
	poly129 buf2;
	poly129 * ptr1 = & buf1;
	poly129 * ptr2 = & buf2;
	poly129 * ptr_tmp;
	///x144_plus_16x[120];
	//*ptr1 = x129;
	//poly129::add( *ptr1 , x129 , x129_ , poly129::gf_t::one() );
	x129.reduce_headterm_add( *ptr1 , x129_ );
	for(int i=15;i>0;i--) {
		ptr1->mul_X();
		ptr1->reduce_term_by( x129 , *ptr2 );
		ptr_tmp=ptr1; ptr1=ptr2; ptr2=ptr_tmp;
	}
	x144_plus_16x[0] = *ptr1;
	for( unsigned i=1;i<120;i++ ) {
		for(int j=16;j>0;j--) {
			ptr1->mul_X();
			ptr1->reduce_term_by( x129 , *ptr2 );
			ptr_tmp=ptr1; ptr1=ptr2; ptr2=ptr_tmp;
		}
		x144_plus_16x[i] = *ptr1;
	}
}


poly129 & poly_pow16_mod_tab::pow16_mod( poly129 & p ) const
{
	poly129 exceed_terms;
	poly129 buf1;
	poly129 buf2;
	poly129 * ptr1 = &buf1;
	poly129 * ptr2 = &buf2;

	p.pow2to(4);
	p.split( exceed_terms , buf1 , 128 );

	for( unsigned i=0; i<exceed_terms.n_terms; i++ ) {
		poly129::add( *ptr2 , *ptr1 , x144_plus_16x[(exceed_terms.deg[i]-144)>>4] , exceed_terms.coef[i] );
		poly129 * tmp = ptr2;
		ptr2 = ptr1;
		ptr1 = tmp;
	}
	p = *ptr1;
	return p;
}


poly129 & poly_pow16_mod_tab::Xpow2to103( poly129 & x2103 ) const
{
	x2103 = x144_plus_16x[119];  // x^2^(7+4)
	for(unsigned i=11;i<103;i+=4) pow16_mod( x2103 );
	return x2103;
}

/////////////////////////////////


poly_pow32_mod_tab::poly_pow32_mod_tab( const poly129 & p )
{
	x129 = p;
	x129.normalize();
	poly129 x129_;
	x129_.n_terms = 1;
	x129_.deg[0] = 129;
	x129_.coef[0] = poly129::gf_t::one();

	poly129 buf1;
	poly129 buf2;
	poly129 * ptr1 = & buf1;
	poly129 * ptr2 = & buf2;
	poly129 * ptr_tmp;
	x129.reduce_headterm_add( *ptr1 , x129_ );
//      poly129 x130_plus_2x[64];
	ptr1->mul_X();
	ptr1->reduce_term_by( x129 , *ptr2 );
	ptr_tmp=ptr1; ptr1=ptr2; ptr2=ptr_tmp;
	x130_plus_2x[0] = *ptr1;
	for( unsigned i=1;i<64;i++) {
		ptr1->mul_X();
		ptr1->reduce_term_by( x129 , *ptr2 );
		ptr_tmp=ptr1; ptr1=ptr2; ptr2=ptr_tmp;
		ptr1->mul_X();
		ptr1->reduce_term_by( x129 , *ptr2 );
		ptr_tmp=ptr1; ptr1=ptr2; ptr2=ptr_tmp;
		x130_plus_2x[i] = *ptr1;
	}
//      poly129 x160_plus_32x[124];
	for( unsigned i=0;i<4;i++)
		x160_plus_32x[i] = x130_plus_2x[15+i*16];

	for( unsigned i=4;i<124;i++ ) {
		for(int j=32;j>0;j--) {
			ptr1->mul_X();
			ptr1->reduce_term_by( x129 , *ptr2 );
			ptr_tmp=ptr1; ptr1=ptr2; ptr2=ptr_tmp;
		}
		x160_plus_32x[i] = *ptr1;
	}
}

poly129 & poly_pow32_mod_tab::squ_mod( poly129 & p ) const
{
	poly129 exceed_terms;
	poly129 buf1;
	poly129 buf2;
	poly129 * ptr1 = &buf1;
	poly129 * ptr2 = &buf2;

	p.squ();
	p.split( exceed_terms , buf1 , 128 );

	for( unsigned i=0; i<exceed_terms.n_terms; i++ ) {
		poly129::add( *ptr2 , *ptr1 , x130_plus_2x[(exceed_terms.deg[i]-130)>>1] , exceed_terms.coef[i] );
		poly129 * tmp = ptr2;
		ptr2 = ptr1;
		ptr1 = tmp;
	}
	p = *ptr1;
	return p;
}

poly129 & poly_pow32_mod_tab::pow32_mod( poly129 & p ) const
{
	poly129 exceed_terms;
	poly129 buf1;
	poly129 buf2;
	poly129 * ptr1 = &buf1;
	poly129 * ptr2 = &buf2;

	p.pow2to(5);
	p.split( exceed_terms , buf1 , 128 );

	for( unsigned i=0; i<exceed_terms.n_terms; i++ ) {
		poly129::add( *ptr2 , *ptr1 , x160_plus_32x[(exceed_terms.deg[i]-160)>>5] , exceed_terms.coef[i] );
		poly129 * tmp = ptr2;
		ptr2 = ptr1;
		ptr1 = tmp;
	}
	p = *ptr1;
	return p;
}

poly129 & poly_pow32_mod_tab::Xpow2to103( poly129 & x2103 ) const
{
	x2103 = x160_plus_32x[123];  // x^2^(7+5)
	squ_mod( x2103 ); // x^2^13
	for(unsigned i=13;i<103;i+=5) pow32_mod( x2103 );
	return x2103;
}






/////////////////////////////////////////////////////////////////////////////////



template <unsigned ext,unsigned max_deg,unsigned pow2 >
struct poly_pow_tab
{
	typedef poly<ext,max_deg> poly_t;
	typedef typename poly_t::gf_t gf_t;

	static const unsigned mod_mask = (1<<pow2)-1;

	unsigned base_deg;
	poly_t x_base; /// zero
	poly_t x_n_table[max_deg];

	poly_t x_deg_table[(1<<pow2)+1];

	poly_pow_tab( const poly_t & nor_p ) {
		x_base = nor_p;
		base_deg = nor_p.deg[0];

		poly_t x_deg;
		x_deg.n_terms = 1;
		x_deg.deg[0] = base_deg;
		x_deg.coef[0] = gf_t::one();

		poly_t buf1;
		poly_t buf2;
		poly_t * ptr1 = & buf1;
		poly_t * ptr2 = & buf2;
		poly_t * tmp;

#ifdef CONFIG_PROFILE
	bm_start(&prepare_poly);
#endif

		poly_t::template _add<true,false,false>( *ptr1 , x_base , x_deg , 0 , gf_t::one() );
#if MAX_DEG < 10
		x_deg_table[0] = *ptr1;
		unsigned i;
		for(i=base_deg+1; i <= base_deg + (1<<pow2); i++ ) {
			ptr1->mul_X();
			ptr1->reduce_term_by( x_base , *ptr2 );
			tmp=ptr1; ptr1=ptr2; ptr2=tmp;
			x_deg_table[i-base_deg] = *ptr1;
			if( 0 == (i&mod_mask) ) x_n_table[i>>pow2] = *ptr1;
		}

		poly_t exceed_terms;
		*ptr1 = x_n_table[i>>pow2];
		for( i = (i&(~mod_mask))+(1<<pow2); i <=((max_deg-1)<<pow2); i += (1<<pow2)) {
			ptr1->mul_X(1<<pow2);
			ptr1->split( exceed_terms , *ptr2 , base_deg-1 );
			tmp = ptr2; ptr2 = ptr1; ptr1 = tmp;
			for( unsigned j=0;j<exceed_terms.n_terms;j++) {
				poly_t::add( *ptr2 , *ptr1 , x_deg_table[exceed_terms.deg[j]-base_deg] , exceed_terms.coef[j] );
				tmp = ptr2; ptr2 = ptr1; ptr1 = tmp;
			}
			x_n_table[i>>pow2] = *ptr1;
		}
#else
		poly_t::template _add<true,false,false>( *ptr1 , x_base , x_deg , 0 , gf_t::one() );
		for(unsigned i=base_deg+1; i <=((max_deg-1)<<pow2) ; i++ ) {
			ptr1->mul_X();
			ptr1->reduce_term_by( x_base , *ptr2 );
			tmp=ptr1; ptr1=ptr2; ptr2=tmp;
			if( 0 == (i&mod_mask) ) x_n_table[i>>pow2] = *ptr1;
		}
#endif
#ifdef CONFIG_PROFILE
	bm_stop(&prepare_poly);
#endif
	}

	poly_t & pow_mod( poly_t & p ) const {
		poly_t exceed_terms;
		poly_t buf1;
		poly_t buf2;
		poly_t * ptr1 = &buf1;
		poly_t * ptr2 = &buf2;
		poly_t * tmp;

		p.pow2to(pow2);
		p.split( exceed_terms , buf1 , base_deg-1 );

		for( unsigned i=0; i<exceed_terms.n_terms; i++ ) {
			poly_t::add( *ptr2 , *ptr1 , x_n_table[exceed_terms.deg[i]>>pow2] , exceed_terms.coef[i] );
			tmp = ptr2; ptr2 = ptr1; ptr1 = tmp;
		}
		p = *ptr1;
		return p;
	}

	void X_2ext( poly_t & p ) const;
};


////////////////////////////



template <>
void poly_pow_tab<103,129,2>::X_2ext( poly_t & p) const
{
#ifdef CONFIG_PROFILE
	BENCHMARK( do_poly , {
#endif

	p = x_n_table[128];  /// X^512 = X^2^9
	for(unsigned i=9;i<103;i+=2) pow_mod( p );

#ifdef CONFIG_PROFILE
	} );
#endif
}

template <>
void poly_pow_tab<103,129,3>::X_2ext( poly_t & p) const
{
#ifdef CONFIG_PROFILE
	BENCHMARK( do_poly , {
#endif

	p = x_n_table[128];  /// X^1024 = X^2^10
	for(unsigned i=10;i<103;i+=3) pow_mod( p );

#ifdef CONFIG_PROFILE
	} );
#endif
}


template <>
void poly_pow_tab<103,129,4>::X_2ext( poly_t & p) const
{
#ifdef CONFIG_PROFILE
	BENCHMARK( do_poly , {
#endif

	p = x_n_table[128];  /// X^(7+4) = X^2^11
	for(unsigned i=11;i<103;i+=4) pow_mod( p );

#ifdef CONFIG_PROFILE
	} );
#endif
}


template <>
void poly_pow_tab<103,129,5>::X_2ext( poly_t & p) const
{
#ifdef CONFIG_PROFILE
	BENCHMARK( do_poly , {
#endif

	p = x_n_table[8];  /// X^(3+5) = X^2^12
	for(unsigned i=8;i<103;i+=5) pow_mod( p );

#ifdef CONFIG_PROFILE
	} );
#endif
}


template <>
void poly_pow_tab<103,129,6>::X_2ext( poly_t & p) const
{
#ifdef CONFIG_PROFILE
	BENCHMARK( do_poly , {
#endif

	p = x_n_table[128];  /// X^(7+6) = X^2^13
	for(unsigned i=13;i<103;i+=6) pow_mod( p );

#ifdef CONFIG_PROFILE
	} );
#endif
}


////////////////////////////////////////////////////////////////


#ifndef USE_OLD

template <>
void poly<103,129>::X_2ext_X( poly_t & x_2ext_x ) const
{
	poly_t x_2ext = *this;
	x_2ext.normalize();

#ifdef CONFIG_PROFILE
	new_phase();
#endif
//	poly_pow_tab<94,17,1> pow_eng( x_2ext );
//	poly_pow_tab<94,17,2> pow_eng( x_2ext );
#ifdef CONFIG_HAS_PCLMULQDQ
	poly_pow_tab<103,129,4> pow_eng( x_2ext );
#else
	poly_pow_tab<103,129,3> pow_eng( x_2ext );
#endif
//	poly_pow_tab<103,129,2> pow_eng( x_2ext );
//	poly_pow_tab<103,129,3> pow_eng( x_2ext );
//	poly_pow_tab<103,129,4> pow_eng( x_2ext );
//	poly_pow_tab<103,129,5> pow_eng( x_2ext );
//	poly_pow_tab<103,129,6> pow_eng( x_2ext );



#ifdef CONFIG_PROFILE
	new_phase();
#endif

	pow_eng.X_2ext( x_2ext );

#ifdef CONFIG_PROFILE
	new_phase();
#endif

	gf_t gf_one = gf_t::one();
	poly_t x;
	x.n_terms = 1;
	x.deg[0] = 1;
	x.coef[0] = gf_one;

	poly_t::add( x_2ext_x , x_2ext , x , gf_one );
}

#endif
