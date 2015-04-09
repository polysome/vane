#include "run_config.h"

#include "rainbow.h"

#include "linear31.h"


static void gen_ov_rowmat( uint32_t r_rowmat[20][20] , const uint8_t ol_rowmat[20][20] , const uint8_t * ov_rowmat 
	, const uint8_t * v_vec , unsigned v_len )
{
	unsigned i,j;
	for(i=0;i<20;i++){
		for(j=0;j<20;j++) {
			r_rowmat[i][j] = ol_rowmat[i][j];
			r_rowmat[i][j] += (uint32_t)vec_dot( ov_rowmat , v_vec , v_len );
			ov_rowmat += v_len;
		}
	}
}




static void sec_pubmap( uint8_t *r , const void *_key , const uint8_t * inp )
{
	const seckey_t * key = (const seckey_t *)_key;
	uint32_t r32[64] = {0};
	uint8_t tmp[64] = {0};
	uint8_t tmp2[40] = {0};
	uint8_t tmp3[20] = {0};

	uint32_t rowmat[20][20];

	vec_assign32( r32 , key->sc , 64 );
	mat_mad32( r32 , &key->s[0][0] , inp , 64 );
	vec_fullreduce32_cvt( tmp , r32 , 64 );

	gen_ov_rowmat( rowmat , key->ol1st_rowmat , & key->ov1st_rowmat[0][0][0] , tmp , 24 );
	rowmat_mul32( tmp3 , &rowmat[0][0] , tmp+24 , 20 );
	eval_q24x20( tmp2 , &(key->vv1st) , tmp );
	vec_add( tmp2 , tmp3 , 20 );

	gen_ov_rowmat( rowmat , key->ol2nd_rowmat , & key->ov2nd_rowmat[0][0][0] , tmp , 44 );
	rowmat_mul32( tmp3 , &rowmat[0][0] , tmp+44 , 20 );
	eval_q44x20( &tmp2[20] , &(key->vv2nd) , tmp );
	vec_add( &tmp2[20] , tmp3 , 20 );

	vec_assign32( r32 , key->tc , 40 );
	mat_mad32( r32 , &key->t[0][0] , tmp2 , 40 );
	vec_fullreduce32_cvt( r , r32 , 40 );
	vec_fullreduce( r , 40 ); /* remove 31 */
}




/*  public functions */

#if !defined(__DEBUG__)
static
#endif
int genkey( uint8_t * pubkey , uint8_t * seckey )
{
	qpoly_64x40_t * pk = (pubkey_t*) pubkey;
	seckey_t * sk = (seckey_t*) seckey;
	uint8_t inv_s[64*64];
	uint8_t inv_t[40*40];
	uint8_t inp64[64]={0};
	uint8_t out40[40]={0};

	vec_rand( sk->sc , 64 );
	vec_rand( (uint8_t *)&sk->vv1st , sizeof(qpoly_24x20_t) );
	vec_rand( &sk->ov1st_rowmat[0][0][0] , 20*20*24 );
	vec_rand( &sk->ol1st_rowmat[0][0] , 20*20 );
	vec_rand( (uint8_t *)&sk->vv2nd , sizeof(qpoly_44x20_t) );
	vec_rand( &sk->ov2nd_rowmat[0][0][0] , 20*20*44 );
	vec_rand( &sk->ol2nd_rowmat[0][0] , 20*20 );
	vec_setzero( sk->tc , 40 );

	mat_rand( sk->s[0] , inv_s , 64 );
	mat_rand( sk->t[0] , inv_t , 40 );

	sec_pubmap( out40 , sk , inp64 );
	vec_negative( sk->tc , out40 , 40 );

	interpolate_64x40( pk , sec_pubmap , (void *)sk );

	vec_assign( sk->s[0] , inv_s , 64*64 );
	vec_assign( sk->t[0] , inv_t , 40*40 );

	vec_negative( sk->sc , sk->sc , 64 );
	vec_negative( sk->tc , sk->tc , 40 );

	return 0;
}

#if !defined(__DEBUG__)
static
#endif
int sign( uint8_t * s , const seckey_t * key , const uint8_t * m )
{
	uint32_t r32[64] = {0};
	uint8_t tmp[64];
	uint8_t tmp2[40];
	uint8_t tmp3[20];
	unsigned i;
	int badluck;

	uint32_t rowmat[20][20];

	vec_assign( tmp2 , m , 40 );
	vec_add( tmp2 , key->tc , 40 );
	mat_mad32( r32 , &key->t[0][0] , tmp2 , 40 );
	vec_fullreduce32_cvt( tmp2 , r32 , 40 );

	for(i=0;i<5;i++) {
		vec_rand(tmp,24);
		eval_q24x20( tmp3 , &(key->vv1st) , tmp );
		vec_negative( tmp3 , tmp3 , 20 );
		vec_add( tmp3 , tmp2 , 20 );
		gen_ov_rowmat( rowmat , key->ol1st_rowmat , & key->ov1st_rowmat[0][0][0] , tmp , 24 );
		badluck = solve_linear20( &tmp[24] , &rowmat[0][0] , tmp3 );
		if( 0 != badluck ) continue;

		eval_q44x20( tmp3 , &(key->vv2nd) , tmp );
		vec_negative( tmp3 , tmp3 , 20 );
		vec_add( tmp3 , &tmp2[20] , 20 );
		gen_ov_rowmat( rowmat , key->ol2nd_rowmat , & key->ov2nd_rowmat[0][0][0] , tmp , 44 );
		badluck = solve_linear20( &tmp[44] , &rowmat[0][0] , tmp3 );
		if( 0 == badluck ) break;
	}
	if( 0 != badluck ) return -1;

	vec_add(tmp, key->sc , 64 );
	vec_setzero((uint8_t *)r32,64*4);
	mat_mad32( r32 , &key->s[0][0] , tmp , 64 );
	vec_fullreduce32_cvt( s , r32 , 64 );

	vec_fullreduce( s , 40 ); /* remove 31 */
	return 0;
}


#if !defined(__DEBUG__)
static
#endif
int verify( const uint8_t * md , const pubkey_t * _key , const uint8_t * s )
{
	const pubkey_t * key = (const pubkey_t *)_key;
	uint8_t r[40];

	eval_q64x40( r , key , s );

	return vec_cmp40(md,r);
}


/*  binary interface  */



int genkey_pack( uint8_t * pubkey , uint8_t *seckey )
{
	int r;
	unsigned i;
	pubkey_t pk;
	uint8_t * bptr = (uint8_t*)&pk;
	r = genkey( bptr , seckey );
	for(i=0;i<sizeof(pubkey_t)/8;i++){
		pack_40b_31x8( pubkey , bptr );
		pubkey += 5;
		bptr += 8;
	}
	return r;
}

int sign_bin( uint8_t * s320b , const seckey_t * key , const uint8_t * md192b )
{
	uint8_t md[40];
	uint8_t s[64];
	int r;
//	vec_dump("sign_bin(): ",md192b, 24 );
	cvt_31x20_bin96(&md[0],(const uint32_t *)&md192b[0]);
	cvt_31x20_bin96(&md[20],(const uint32_t *)&md192b[12]);
//	vec_dump("cvt: ",md,40);
	r = sign( s , key , md );
//	vec_dump("sign->:",s,64);
	pack_40b_31x8( s320b , s );
	pack_40b_31x8( &s320b[5] , &s[8] );
	pack_40b_31x8( &s320b[10] , &s[16] );
	pack_40b_31x8( &s320b[15] , &s[24] );
	pack_40b_31x8( &s320b[20] , &s[32] );
	pack_40b_31x8( &s320b[25] , &s[40] );
	pack_40b_31x8( &s320b[30] , &s[48] );
	pack_40b_31x8( &s320b[35] , &s[56] );
//	vec_dump("pack->",s320b,40);
	return r;
}


inline static void vec_mad32( uint32_t * accu_r , const uint8_t * vec , unsigned c , unsigned len )
{
	unsigned i;
	for(i=0;i<len;i++) accu_r[i] += c*((unsigned)vec[i]);
}
inline static void vec_mul32( uint32_t * r , const uint8_t * vec , unsigned c , unsigned len )
{
	unsigned i;
	for(i=0;i<len;i++) r[i] = c*((unsigned)vec[i]);
}


int verify_bin( const uint8_t * md192b , const uint8_t * key , const uint8_t * s320b )
{
	uint8_t s[64];
	uint32_t accu_r[40] = {0};
	uint32_t tmp[64];
	unsigned i,j;
	uint8_t partial_key[40];

//	vec_dump("verify():",s320b,40);
	unpack_31x8_40b( s , s320b );
	unpack_31x8_40b( &s[8] , &s320b[5] );
	unpack_31x8_40b( &s[16] , &s320b[10] );
	unpack_31x8_40b( &s[24] , &s320b[15] );
	unpack_31x8_40b( &s[32] , &s320b[20] );
	unpack_31x8_40b( &s[40] , &s320b[25] );
	unpack_31x8_40b( &s[48] , &s320b[30] );
	unpack_31x8_40b( &s[56] , &s320b[35] );

//	vec_dump("s: ",s,64);

        for(i=0;i<64;i++){
		unpack_31x8_40b(&partial_key[0],&key[0]);
		unpack_31x8_40b(&partial_key[8],&key[5]);
		unpack_31x8_40b(&partial_key[16],&key[10]);
		unpack_31x8_40b(&partial_key[24],&key[15]);
		unpack_31x8_40b(&partial_key[32],&key[20]);
		key += 25;
		vec_mad32(accu_r,partial_key,s[i],40);
	}
        for(i=0;i<64;i++){
		vec_mul32( tmp , &s[i] , s[i] , 64-i );
		for(j=0;j<64-i;j++) {
			unpack_31x8_40b(&partial_key[0],&key[0]);
			unpack_31x8_40b(&partial_key[8],&key[5]);
			unpack_31x8_40b(&partial_key[16],&key[10]);
			unpack_31x8_40b(&partial_key[24],&key[15]);
			unpack_31x8_40b(&partial_key[32],&key[20]);
			key += 25;
			vec_mad32(accu_r,partial_key,tmp[j],40);
		}
	}
	vec_fullreduce32_cvt( s , accu_r , 40 );
	vec_fullreduce( s , 40 );

//	vec_dump("md:",s,40);

	cvt_bin96_31x20( (uint32_t *) &partial_key[0] , s );
	cvt_bin96_31x20( (uint32_t *) &partial_key[12] , &s[20] );

//	vec_dump("cvt:",partial_key,24);

	j = 0;
	for( i=0;i<6;i++) j |= ((uint32_t*)partial_key)[i]^((const uint32_t *)md192b)[i];

	return (0==j)?0:-1;
}




