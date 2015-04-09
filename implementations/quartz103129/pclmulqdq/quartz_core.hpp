
#ifndef _QUARTZ_CORE_HPP_
#define _QUARTZ_CORE_HPP_

#include "mpkc.hpp"
#include "gf2ext_poly.hpp"


/// pub key  MAT<TERMS_QUAD_POLY<N>,M>


template <unsigned ext,unsigned max_d,unsigned minus,unsigned vinegar>
struct quartz_key{
	typedef GF2EXT<ext> gf_t;
	typedef poly<ext,max_d> poly_t;

	static const unsigned t_size = ext;
	static const unsigned s_size = ext+vinegar;
	static const unsigned roll_size = minus + vinegar;

	void gen_univar_poly( poly_t & p , const VEC<vinegar>& v ) const;

};

template <unsigned ext,unsigned minus,unsigned vinegar>
struct quartz_key<ext,129,minus,vinegar>{
	static const unsigned max_d = 129;
	//static const unsigned minus = 3;
	//static const unsigned vinegar = 4;

	typedef GF2EXT<ext> gf_t;
	typedef poly<ext,max_d> poly_t;

	static const unsigned t_size = ext;
	static const unsigned s_size = ext+vinegar;
	static const unsigned roll_size = minus + vinegar;

	typedef MAT<t_size,t_size> mat_t_t;
	typedef VEC<t_size> vec_t_t;
	typedef MAT<s_size,s_size> mat_s_t;
	typedef VEC<s_size> vec_s_t;

	MAT<t_size,t_size> t_mat;
	VEC<t_size> t_vec;

	MAT<s_size,s_size> s_mat;
	VEC<s_size> s_vec;

	gf_t q_a_2[1];
	gf_t q_a_4[2];
	gf_t q_a_8[3];
	gf_t q_a_16[4];
	gf_t q_a_32[5];
	gf_t q_a_64[6];
	gf_t q_a_128[1];

	gf_t l_b_1[vinegar];
	gf_t l_b_2[vinegar];
	gf_t l_b_4[vinegar];
	gf_t l_b_8[vinegar];
	gf_t l_b_16[vinegar];
	gf_t l_b_32[vinegar];
	gf_t l_b_64[vinegar];
	gf_t l_b_128[vinegar];

	gf_t c_r[vinegar*vinegar];

	static unsigned num_byte() { return mat_t_t::num_byte()+vec_t_t::num_byte()+mat_s_t::num_byte()+vec_s_t::num_byte()+gf_t::num_byte()*(22+ (8+vinegar)*vinegar);}
	void dump( uint8_t * m ) const {
		t_mat.dump(m); m += mat_t_t::num_byte();
		t_vec.dump(m); m += vec_t_t::num_byte();
		s_mat.dump(m); m += mat_s_t::num_byte();
		s_vec.dump(m); m += vec_s_t::num_byte();
		q_a_2[0].dump(m); m += gf_t::num_byte();
		for(unsigned i=0;i<2;i++) { q_a_4[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<3;i++) { q_a_8[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<4;i++) { q_a_16[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<5;i++) { q_a_32[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<6;i++) { q_a_64[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<1;i++) { q_a_128[i].dump(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar;i++) { l_b_1[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_2[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_4[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_8[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_16[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_32[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_64[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_128[i].dump(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar*vinegar;i++) { c_r[i].dump(m); m += gf_t::num_byte(); }
	}
	void set( const uint8_t * m ) {
		t_mat.set(m); m += mat_t_t::num_byte();
		t_vec = vec_t_t(m); m += vec_t_t::num_byte();
		s_mat.set(m); m += mat_s_t::num_byte();
		s_vec = vec_s_t(m); m += vec_s_t::num_byte();

		q_a_2[0] = gf_t(m); m += gf_t::num_byte();
		for(unsigned i=0;i<2;i++) { q_a_4[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<3;i++) { q_a_8[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<4;i++) { q_a_16[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<5;i++) { q_a_32[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<6;i++) { q_a_64[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<1;i++) { q_a_128[i] = gf_t(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar;i++) { l_b_1[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_2[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_4[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_8[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_16[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_32[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_64[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_128[i] = gf_t(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar*vinegar;i++) { c_r[i] = gf_t(m); m += gf_t::num_byte(); }
	}

	void set_central_random() {
		q_a_2[0] = gf_t::rand();
		for(int i=0;i<2;i++) q_a_4[i] = gf_t::rand();
		for(int i=0;i<3;i++) q_a_8[i] = gf_t::rand();
		for(int i=0;i<4;i++) q_a_16[i] = gf_t::rand();
		for(int i=0;i<5;i++) q_a_32[i] = gf_t::rand();
		for(int i=0;i<6;i++) q_a_64[i] = gf_t::rand();
		q_a_128[0] = gf_t::rand();

		for(unsigned i=0;i<vinegar;i++) l_b_1[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_1[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_2[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_4[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_8[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_16[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_32[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_64[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_128[i] = gf_t::rand();

		for(unsigned i=0;i<vinegar*vinegar;i++) c_r[i] = gf_t::rand();
	}

	void get_central_poly( poly_t & r_poly , uint32_t v ) const {
		gf_t cc;
		r_poly.set_zero();

		r_poly.append_term(q_a_128[0],129);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_128[i]:gf_t::zero();
		r_poly.append_term(cc,128);

		r_poly.append_term(q_a_64[5],96);
		r_poly.append_term(q_a_64[4],80);
		r_poly.append_term(q_a_64[3],72);
		r_poly.append_term(q_a_64[2],68);
		r_poly.append_term(q_a_64[1],66);
		r_poly.append_term(q_a_64[0],65);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_64[i]:gf_t::zero();
		r_poly.append_term(cc,64);

		r_poly.append_term(q_a_32[4],48);
		r_poly.append_term(q_a_32[3],40);
		r_poly.append_term(q_a_32[2],36);
		r_poly.append_term(q_a_32[1],34);
		r_poly.append_term(q_a_32[0],33);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_32[i]:gf_t::zero();
		r_poly.append_term(cc,32);

		r_poly.append_term(q_a_16[3],24);
		r_poly.append_term(q_a_16[2],20);
		r_poly.append_term(q_a_16[1],18);
		r_poly.append_term(q_a_16[0],17);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_16[i]:gf_t::zero();
		r_poly.append_term(cc,16);

		r_poly.append_term(q_a_8[2],12);
		r_poly.append_term(q_a_8[1],10);
		r_poly.append_term(q_a_8[0],9);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_8[i]:gf_t::zero();
		r_poly.append_term(cc,8);

		r_poly.append_term(q_a_4[1],6);
		r_poly.append_term(q_a_4[0],5);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_4[i]:gf_t::zero();
		r_poly.append_term(cc,4);

		r_poly.append_term(q_a_2[0],3);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_2[i]:gf_t::zero();
		r_poly.append_term(cc,2);

		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_1[i]:gf_t::zero();
		r_poly.append_term(cc,1);

		cc.set_zero();
		for(unsigned i=0;i<vinegar;i++) {
			bool add = (1<<i)&v;
			for(unsigned j=0;j<vinegar;j++) {
				bool add2 = (1<<j) & v;
				cc ^= (add & add2) ? c_r[i*vinegar+j]:gf_t::zero();
			}
		}
		r_poly.append_term(cc,0);
	}

};


template <unsigned ext,unsigned minus,unsigned vinegar>
struct quartz_key<ext,17,minus,vinegar>{
	//static const unsigned ext = 94;
	static const unsigned max_d = 17;
	//static const unsigned minus = 4;
	//static const unsigned vinegar = 4;

	typedef GF2EXT<ext> gf_t;
	typedef poly<ext,max_d> poly_t;

	static const unsigned t_size = ext;
	static const unsigned s_size = ext+vinegar;
	static const unsigned roll_size = minus + vinegar;

	typedef MAT<t_size,t_size> mat_t_t;
	typedef VEC<t_size> vec_t_t;
	typedef MAT<s_size,s_size> mat_s_t;
	typedef VEC<s_size> vec_s_t;

	MAT<t_size,t_size> t_mat;
	VEC<t_size> t_vec;

	MAT<s_size,s_size> s_mat;
	VEC<s_size> s_vec;

	gf_t q_a_2[1];
	gf_t q_a_4[2];
	gf_t q_a_8[3];
	gf_t q_a_16[1];

	gf_t l_b_1[vinegar];
	gf_t l_b_2[vinegar];
	gf_t l_b_4[vinegar];
	gf_t l_b_8[vinegar];
	gf_t l_b_16[vinegar];

	gf_t c_r[vinegar*vinegar];

	static unsigned num_byte() { return mat_t_t::num_byte()+vec_t_t::num_byte()+mat_s_t::num_byte()+vec_s_t::num_byte()+gf_t::num_byte()*(7+ (5+vinegar)*vinegar);}
	void dump( uint8_t * m ) const {
		t_mat.dump(m); m += mat_t_t::num_byte();
		t_vec.dump(m); m += vec_t_t::num_byte();
		s_mat.dump(m); m += mat_s_t::num_byte();
		s_vec.dump(m); m += vec_s_t::num_byte();
		q_a_2[0].dump(m); m += gf_t::num_byte();
		for(unsigned i=0;i<2;i++) { q_a_4[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<3;i++) { q_a_8[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<1;i++) { q_a_16[i].dump(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar;i++) { l_b_1[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_2[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_4[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_8[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_16[i].dump(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar*vinegar;i++) { c_r[i].dump(m); m += gf_t::num_byte(); }
	}
	void set( const uint8_t * m ) {
		t_mat.set(m); m += mat_t_t::num_byte();
		t_vec = vec_t_t(m); m += vec_t_t::num_byte();
		s_mat.set(m); m += mat_s_t::num_byte();
		s_vec = vec_s_t(m); m += vec_s_t::num_byte();

		q_a_2[0] = gf_t(m); m += gf_t::num_byte();
		for(unsigned i=0;i<2;i++) { q_a_4[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<3;i++) { q_a_8[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<1;i++) { q_a_16[i] = gf_t(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar;i++) { l_b_1[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_2[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_4[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_8[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_16[i] = gf_t(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar*vinegar;i++) { c_r[i] = gf_t(m); m += gf_t::num_byte(); }
	}

	void set_central_random() {
		q_a_2[0] = gf_t::rand();
		for(int i=0;i<2;i++) q_a_4[i] = gf_t::rand();
		for(int i=0;i<3;i++) q_a_8[i] = gf_t::rand();
		q_a_16[0] = gf_t::rand();

		for(unsigned i=0;i<vinegar;i++) l_b_1[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_2[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_4[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_8[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_16[i] = gf_t::rand();

		for(unsigned i=0;i<vinegar*vinegar;i++) c_r[i] = gf_t::rand();
	}

	void get_central_poly( poly_t & r_poly , uint32_t v ) const {
		gf_t cc;
		r_poly.set_zero();

		r_poly.append_term(q_a_16[0],17);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_16[i]:gf_t::zero();
		r_poly.append_term(cc,16);

		r_poly.append_term(q_a_8[2],12);
		r_poly.append_term(q_a_8[1],10);
		r_poly.append_term(q_a_8[0],9);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_8[i]:gf_t::zero();
		r_poly.append_term(cc,8);

		r_poly.append_term(q_a_4[1],6);
		r_poly.append_term(q_a_4[0],5);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_4[i]:gf_t::zero();
		r_poly.append_term(cc,4);

		r_poly.append_term(q_a_2[0],3);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_2[i]:gf_t::zero();
		r_poly.append_term(cc,2);

		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_1[i]:gf_t::zero();
		r_poly.append_term(cc,1);

		cc.set_zero();
		for(unsigned i=0;i<vinegar;i++) {
			bool add = (1<<i)&v;
			for(unsigned j=0;j<vinegar;j++) {
				bool add2 = (1<<j) & v;
				cc ^= (add & add2) ? c_r[i*vinegar+j]:gf_t::zero();
			}
		}
		r_poly.append_term(cc,0);
	}

};


template <unsigned ext,unsigned minus,unsigned vinegar>
struct quartz_key<ext,9,minus,vinegar>{
	//static const unsigned ext = 95;
	static const unsigned max_d = 9;
	//static const unsigned minus = 5;
	//static const unsigned vinegar = 5;

	typedef GF2EXT<ext> gf_t;
	typedef poly<ext,max_d> poly_t;

	static const unsigned t_size = ext;
	static const unsigned s_size = ext+vinegar;
	static const unsigned roll_size = minus + vinegar;

	typedef MAT<t_size,t_size> mat_t_t;
	typedef VEC<t_size> vec_t_t;
	typedef MAT<s_size,s_size> mat_s_t;
	typedef VEC<s_size> vec_s_t;

	MAT<t_size,t_size> t_mat;
	VEC<t_size> t_vec;

	MAT<s_size,s_size> s_mat;
	VEC<s_size> s_vec;

	gf_t q_a_2[1];
	gf_t q_a_4[2];
	gf_t q_a_8[1];

	gf_t l_b_1[vinegar];
	gf_t l_b_2[vinegar];
	gf_t l_b_4[vinegar];
	gf_t l_b_8[vinegar];

	gf_t c_r[vinegar*vinegar];

	static unsigned num_byte() { return mat_t_t::num_byte()+vec_t_t::num_byte()+mat_s_t::num_byte()+vec_s_t::num_byte()+gf_t::num_byte()*(4+ (4+vinegar)*vinegar);}
	void dump( uint8_t * m ) const {
		t_mat.dump(m); m += mat_t_t::num_byte();
		t_vec.dump(m); m += vec_t_t::num_byte();
		s_mat.dump(m); m += mat_s_t::num_byte();
		s_vec.dump(m); m += vec_s_t::num_byte();
		q_a_2[0].dump(m); m += gf_t::num_byte();
		for(unsigned i=0;i<2;i++) { q_a_4[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<1;i++) { q_a_8[i].dump(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar;i++) { l_b_1[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_2[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_4[i].dump(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_8[i].dump(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar*vinegar;i++) { c_r[i].dump(m); m += gf_t::num_byte(); }
	}
	void set( const uint8_t * m ) {
		t_mat.set(m); m += mat_t_t::num_byte();
		t_vec = vec_t_t(m); m += vec_t_t::num_byte();
		s_mat.set(m); m += mat_s_t::num_byte();
		s_vec = vec_s_t(m); m += vec_s_t::num_byte();

		q_a_2[0] = gf_t(m); m += gf_t::num_byte();
		for(unsigned i=0;i<2;i++) { q_a_4[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<1;i++) { q_a_8[i] = gf_t(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar;i++) { l_b_1[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_2[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_4[i] = gf_t(m); m += gf_t::num_byte(); }
		for(unsigned i=0;i<vinegar;i++) { l_b_8[i] = gf_t(m); m += gf_t::num_byte(); }

		for(unsigned i=0;i<vinegar*vinegar;i++) { c_r[i] = gf_t(m); m += gf_t::num_byte(); }
	}

	void set_central_random() {
		q_a_2[0] = gf_t::rand();
		for(int i=0;i<2;i++) q_a_4[i] = gf_t::rand();
		q_a_8[0] = gf_t::rand();

		for(unsigned i=0;i<vinegar;i++) l_b_1[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_2[i] = gf_t::rand();
		for(unsigned  i=0;i<vinegar;i++) l_b_4[i] = gf_t::rand();
		for(unsigned i=0;i<vinegar;i++) l_b_8[i] = gf_t::rand();

		for(unsigned i=0;i<vinegar*vinegar;i++) c_r[i] = gf_t::rand();
	}

	void get_central_poly( poly_t & r_poly , uint32_t v ) const {
		gf_t cc;
		r_poly.set_zero();

		r_poly.append_term(q_a_8[0],9);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_8[i]:gf_t::zero();
		r_poly.append_term(cc,8);

		r_poly.append_term(q_a_4[1],6);
		r_poly.append_term(q_a_4[0],5);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_4[i]:gf_t::zero();
		r_poly.append_term(cc,4);

		r_poly.append_term(q_a_2[0],3);
		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_2[i]:gf_t::zero();
		r_poly.append_term(cc,2);

		cc.set_zero(); for(unsigned i=0;i<vinegar;i++) cc ^= ((1<<i)&v)?l_b_1[i]:gf_t::zero();
		r_poly.append_term(cc,1);

		cc.set_zero();
		for(unsigned i=0;i<vinegar;i++) {
			bool add = (1<<i)&v;
			for(unsigned j=0;j<vinegar;j++) {
				bool add2 = (1<<j) & v;
				cc ^= (add & add2) ? c_r[i*vinegar+j]:gf_t::zero();
			}
		}
		r_poly.append_term(cc,0);
	}

};




template <unsigned ext,unsigned max_d,unsigned minus,unsigned vinegar>
void quartz_pubmap_seckey( VEC<ext>& z , const quartz_key<ext,max_d,minus,vinegar> & sk , const VEC<ext+vinegar> & w )
//void quartz_pubmap_seckey( VEC<ext-minus>& z , const quartz_key<ext,max_d,minus,vinegar> & sk , const VEC<ext+vinegar> & w )
{
	typedef typename quartz_key<ext,max_d,minus,vinegar>::poly_t poly_t;
	typedef typename quartz_key<ext,max_d,minus,vinegar>::gf_t gf_t;

	VEC<ext+vinegar> qq = sk.s_mat.prod(w);
	qq ^= sk.s_vec;

	unsigned v=0;
	for(unsigned i=0;i<vinegar;i++) v |= ((unsigned)qq.get_ele(ext+i))<<i;

	poly_t cp;
	sk.get_central_poly(cp,v);

	gf_t x( (const uint8_t *)&qq );
	gf_t yy = cp.eval(x);

	VEC<ext> y( (const uint8_t *) &yy );
	VEC<ext> zz = sk.t_mat.prod(y);
	zz ^= sk.t_vec;

	z = zz;

}

template <unsigned ext,unsigned max_d,unsigned minus,unsigned vinegar>
void quartz_wrapper( void * z , const void * key , const void * w )
{
	VEC<ext> _z;

	quartz_pubmap_seckey<ext,max_d,minus,vinegar>( _z, *(const quartz_key<ext,max_d,minus,vinegar>*)key , *(const VEC<ext+vinegar>*)w );

	*(VEC<ext-minus>*) z = _z. template concate<ext-minus>();
//	*(VEC<ext>*) z = _z;
}

template <unsigned ext,unsigned max_d,unsigned minus,unsigned vinegar>
void quartz_gen_key( MAT<TERMS_QUAD_POLY(ext+vinegar),ext-minus>& pk , quartz_key<ext,max_d,minus,vinegar> & sk )
//void quartz_gen_key( MAT<TERMS_QUAD_POLY(ext+vinegar),ext>& pk , quartz_key<ext,max_d,minus,vinegar> & sk )
{
	typedef quartz_key<ext,max_d,minus,vinegar> sk_t;

	MAT<sk_t::t_size,sk_t::t_size> inv_t_mat;
	MAT<sk_t::s_size,sk_t::s_size> inv_s_mat;

	MAT<sk_t::t_size,sk_t::t_size>::rand_inv(inv_t_mat,sk.t_mat);
	MAT<sk_t::s_size,sk_t::s_size>::rand_inv(inv_s_mat,sk.s_mat);

	sk.s_vec = VEC<sk_t::s_size>::rand();

	sk.set_central_random();

	VEC<sk_t::s_size> inp;
	VEC<sk_t::t_size> r;
	inp.set_zero();
	quartz_pubmap_seckey(r,sk,inp);
	sk.t_vec = r;

	interpolate<ext+vinegar,ext-minus>( pk , quartz_wrapper<ext,max_d,minus,vinegar> , &sk );
//	interpolate<ext+vinegar,ext>( pk , quartz_wrapper<ext,max_d,minus,vinegar> , &sk );


	sk.t_mat = inv_t_mat;
	sk.s_mat = inv_s_mat;
}



#include "crypto_hash_sha256.h"
#include <string.h>

#define _LEN_SHA256_ 32

template <unsigned ext,unsigned max_d,unsigned minus,unsigned vinegar>
void quartz_sec_map( VEC<ext+vinegar>& w , const quartz_key<ext,max_d,minus,vinegar> & sk , const VEC<ext-minus> & z , const uint8_t seed[_LEN_SHA256_])
{
	typedef quartz_key<ext,max_d,minus,vinegar> sk_t;
	typedef typename sk_t::gf_t gf_t;
	typedef typename sk_t::poly_t poly_t;

	static uint32_t sha256_digest[_LEN_SHA256_/4];
	memcpy( sha256_digest , seed , _LEN_SHA256_ );

	uint32_t vin = 0;
	const unsigned mask_vin = (1<<vinegar)-1;
	uint32_t min = 0;
	const unsigned mask_min = (1<<minus)-1;

	if( z.is_zero() ) { w.set_zero(); return; }

	VEC<ext> _z;
	VEC<ext> y;

	gf_t gfy;
	gf_t gfx;
	poly_t s_poly;

	bool uni_root = false;
	for(int i=0;i<(2<<vinegar);i++) {
/// roll minus and vinegar /// set minus part constant ???
///		RAND_bytes( (unsigned char*) &vin , 4 );
///		SHA256( (const uint8_t *)sha256_digest , SHA256_DIGEST_LENGTH , (uint8_t *)sha256_digest );
		crypto_hash_sha256( (uint8_t *)sha256_digest ,  (const uint8_t *)sha256_digest , _LEN_SHA256_ );
		vin = sha256_digest[0];

		min = vin & mask_min;
		vin >>= minus; vin &= mask_vin;

		_z = z. template concate<ext>( min );
		_z ^= sk.t_vec;
		y = sk.t_mat.prod(_z);

		gfy = gf_t( (const uint8_t *) &y );
		sk.get_central_poly( s_poly , vin );
		s_poly += gfy;
		uni_root = s_poly.find_unique_root( gfx );
		if(uni_root) break;
	};

	VEC<ext> _x( (const uint8_t *)&gfx );
	VEC<ext+vinegar> x = _x. template concate<ext+vinegar>( vin );
	x ^= sk.s_vec;
	w = sk.s_mat.prod( x );
	if(!uni_root) w.set_zero();
}



#endif /// _QUARTZ_CORE_HPP_
