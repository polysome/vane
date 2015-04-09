
#include "sizes.h"

#include <openssl/rand.h>

#include "quartz.h"

#include "quartz.hpp"

//extern "C"
int keypair( unsigned char sk[SECRETKEY_BYTES] , unsigned long long * sklen , unsigned char pk[PUBLICKEY_BYTES] , unsigned long long * pklen )
{
	quartz_pub_key_t pkey;
	quartz_sec_key_t skey;

	quartz_gen_key(pkey,skey);

	pkey.dump(pk);
	skey.dump(sk);

	*sklen = quartz_sec_key_t::num_byte();
	*pklen = quartz_pub_key_t::num_byte();

	return 0;
}


//extern "C"
int signatureofshorthash( unsigned char sm[SIGNATURE_BYTES],unsigned long long *smlen,
	const unsigned char m[SHORTHASH_BYTES],const unsigned long long mlen,
	const unsigned char sk[SECRETKEY_BYTES],const unsigned long long sklen )
{
	if( sklen != SECRETKEY_BYTES ) return -11;
	if (mlen != SHORTHASH_BYTES) return -12;

	quartz_sec_key_t skey;
	skey.set( sk );

	vec_sign_t signature;
	quartz_sign<REPEAT>( signature , m , skey );

	signature.dump(sm);
	*smlen = vec_sign_t::num_byte();

	return 0;
}


//extern "C"
int verification( const unsigned char m[SHORTHASH_BYTES],const unsigned long long mlen,
	const unsigned char sm[SIGNATURE_BYTES],const unsigned long long smlen,
	const unsigned char pk[PUBLICKEY_BYTES],const unsigned long long pklen )
{
	if (smlen != SIGNATURE_BYTES) return -101;
	if (mlen != SHORTHASH_BYTES) return -102;
	if( pklen != PUBLICKEY_BYTES ) return -103;

	quartz_pub_key_t pkey;
	pkey.set( pk );

	vec_sign_t signature(sm);

	return quartz_verify<REPEAT>( m , sm , pkey );
}
