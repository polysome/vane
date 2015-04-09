#ifndef _QUARTZ_CORE_H_
#define _QUARTZ_CORE_H_



#include "config.h"
#include "quartz_core.hpp"



typedef struct quartz_key<CORE_SIZE,MAX_DEG,MINUS,VINEGAR> quartz_sec_key_t;

typedef MAT<TERMS_QUAD_POLY(N),M> quartz_pub_key_t;

typedef VEC<N> vec_n_t;

typedef VEC<M> vec_m_t;



//#define SIZE_BYTE_SEC_KEY (sizeof(quartz_sec_key_t))

//#define SIZE_BYTE_PUBKEY (sizeof(quartz_pub_key_t))


#endif
