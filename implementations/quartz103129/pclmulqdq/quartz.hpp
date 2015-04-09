
#ifndef _QUARTZ_HPP_
#define _QUARTZ_HPP_


#include "quartz_core.h"


template <unsigned width, unsigned times> inline
void split_hash( VEC<width> hn[] , const uint8_t *sha256_dig )
{
	typedef VEC<width> vec_t;
	const uint64_t * dig1 = (const uint64_t *)sha256_dig;
	uint64_t dig2[8] __attribute__((aligned(32))) = {0}; /// 32 byte = 256 bit
	crypto_hash_sha256( (uint8_t *)dig2 , (const uint8_t *)sha256_dig , 32 );

	hn[0] = vec_t(sha256_dig);
	hn[1] = vec_t(&sha256_dig[16]);
	hn[2] = vec_t((const uint8_t*)dig2);
	if(3<times) hn[3] = vec_t(&((const uint8_t*)dig2)[16]);
}


#include "quartz_core.h"



typedef VEC<M+REPEAT*(MINUS+VINEGAR)> vec_sign_t;




template <unsigned times>
int quartz_verify( const unsigned char * hash256, const vec_sign_t & sm, const quartz_pub_key_t & pk )
{
	vec_m_t hn[times];
	split_hash<M,times>( hn , hash256 );

	vec_n_t nn = sm. template concate<N>();
	vec_m_t accu_check;
	mpkc_pub_map( accu_check , pk , nn );

	uint64_t tail = sm.template tail<(times-1)*(MINUS+VINEGAR)>();
	for(unsigned i=times-1;i>0;i--){
		accu_check ^= hn[i];
		nn = accu_check. template concate<N>( tail );
		mpkc_pub_map( accu_check , pk , nn );
		tail >>= (MINUS+VINEGAR);
	}
	if( accu_check == hn[0] ) return 0;
	return -1;
}



template <unsigned times>
int quartz_sign( vec_sign_t & sm , const unsigned char * hash256, const quartz_sec_key_t & sk )
{
	vec_m_t hn[times];
	split_hash<M,times>( hn , hash256 );

	uint8_t rand_seed[32];
	uint64_t tail=0;

	vec_n_t nn;
	vec_m_t accu_mm;

	accu_mm.set_zero();

	for(unsigned i=0;i<times;i++) {
		accu_mm ^= hn[i];
		memset( rand_seed , 0 , 32 );
		((vec_m_t &)rand_seed) ^= accu_mm;

		quartz_sec_map( nn , sk , accu_mm , rand_seed );

		uint64_t tmp = nn. template tail<MINUS+VINEGAR>();
		tail = (tail<<(MINUS+VINEGAR)) | tmp;
		accu_mm = nn. template concate<M>();
	}

	sm = accu_mm. template concate<M+times*(MINUS+VINEGAR)>(tail);

	return 0;
}



#endif /// _QUARTZ_HPP_

