#if 0

#include "quartz.hpp"
#include <string.h>
#include "benchmark.h"

const unsigned ext = CORE_SIZE;
const unsigned max_d = MAX_DEG;
const unsigned minus = MINUS;
const unsigned vinegar = VINEGAR;

const unsigned num_n = ext+vinegar;
const unsigned num_m = ext-minus;

#define TEST_RUN 50



int main()
{
	union{
	uint64_t s64[4];
	uint8_t sha_seed[_LEN_SHA256_];
	};
	s64[0] = 0;
	s64[1] = 0;
	s64[2] = 0;
	s64[3] = 0;

	quartz_pub_key_t pk;
	quartz_sec_key_t sk;

	quartz_gen_key(pk,sk);


	printf("test:  quartz<%d,%d,%d,%d>\n\n", ext , max_d , minus , vinegar );
	printf("repeat: %d\n",REPEAT);
	printf("N: %d , M: %d\n", num_n , num_m );
	printf("signature size: %d bit\n", ext-minus+(minus+vinegar)*REPEAT);


	vec_sign_t sm;
	sha_seed[0] = 1;
	sha_seed[1] = 2;

	printf("message: ");
	printf("%llx,%llx,%llx,%llx\n",s64[3],s64[2],s64[1],s64[0]);

	int sr =  quartz_sign<REPEAT>( sm , sha_seed , sk );

	printf("signature: ");
	sm.fdump(stdout);
	printf("\n");

	int vr = quartz_verify<REPEAT>( sha_seed , sm , pk );

	printf("verify: %s\n\n", (vr)?"No":"Yes" );


	printf("\nRunning %d tests: ", TEST_RUN);

	benchmark bm_sign,bm_verify;
	bm_init( & bm_sign );
	bm_init( & bm_verify );


	int count = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		RAND_bytes( sha_seed , 32 );
BENCHMARK( bm_sign , {
		quartz_sign<REPEAT>( sm , sha_seed , sk );
});
		int tt;
BENCHMARK( bm_verify , {
		tt = quartz_verify<REPEAT>( sha_seed , sm , pk );
});
		if( 0 == tt ) count++;
	}

	printf("%d / %d passed\n", count , TEST_RUN );
	printf("\n\n");

	char msg[256];
	bm_dump( msg , 256 , & bm_sign );
	printf("sign():\n%s\n",msg);
	bm_dump( msg , 256 , & bm_verify );
	printf("verify():\n%s\n",msg);
	printf("\n\n");

	return 0;
}


#endif
