#ifndef _LINEAR_31_H_
#define _LINEAR_31_H_


#ifdef __cplusplus
extern "C" {
#endif

#include "run_config.h"


#if defined(__DEBUG__)
#include <stdio.h>
#endif


typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned uint32_t;
typedef unsigned long long uint64_t;

/* GF(31) library */

static inline uint8_t fast_mod31_32( unsigned a )
{
	unsigned r = a & 0x7fff;
	r += (a>>15);
	r = (r>>10) + ( r&0x3ff);
	r = (r>>5) + ( r&0x1f);
	return r;
}

static inline uint8_t full_mod31_32( unsigned a )
{
#if 0
	uint8_t r = fast_mod31_32(a); /* 6 bit */
	r = (r>>5) + (r&0x1f);  /* 0 - 32 */
	r = (r>>5) + (r&0x1f);  /* 0 - 31 */
#else
	uint32_t t = ((0x84210842ULL * a)>>36);
	uint8_t r = a + t - (t<<5);
#endif
#if defined(__DEBUG__)
	if(r>=32) printf("full_mod31_32()->%d\n",r);
#endif
	return r;
}

static inline uint32_t fast_mod31_16( uint32_t a )
{
	return (a&0x1f) + (a>>5);
}

static inline uint8_t full_mod31_16( uint32_t a )
{
	unsigned r = (a&0x3ff) + (a>>10);
	r = (r>>5) + (r&0x1f);  /* 0 - 32 */
	r = (r>>5) + (r&0x1f);  /* 0 - 31 */
#if defined(__DEBUG__)
	if(r>=32) printf("full_mod31_16()->%d\n",r);
#endif
	return r;
}


static inline uint8_t full_mod31( uint8_t a )
{
	uint8_t r=a;
	r = (r>>5) + (r&0x1f);  /* 0 - 32 */
	r = (r>>5) + (r&0x1f);  /* 0 - 31 */
#if defined(__DEBUG__)
	if(r>=32) printf("full_mod31()->%d\n",r);
#endif
	return r;
}

static inline uint8_t negative_31( uint8_t a )
{
#if defined(__DEBUG__)
	if(a>=32) printf("negative(): %d->??\n",a);
#endif
	return a^0x1f;
}

static inline uint8_t inv_31( uint8_t a )
{
	static const uint8_t _inv[32] = { 0, 1, 16, 21, 8, 25, 26, 9, 4, 7, 28, 17, 13, 12, 20, 29,
		2, 11, 19, 18, 14, 3, 24, 27, 22, 5, 6, 23, 10, 15, 30, 0 };
#if defined(__DEBUG__)
	if(a>=32) printf("inv_31(): %d->??\n",a);
#endif
	return _inv[a&0x1f];
}

static inline uint8_t remove_31( uint8_t a )
{
	int qq = a;
#if defined(__DEBUG__)
	if(a>=32) printf("remove_31(): %d->??\n",a);
#endif
	qq -= 31;  /* 31->0 , 0-31 -> -??? */
	return (a&(qq>>5));
}



/* binary <--> 31 conversion */

static inline void div31(unsigned * q , unsigned *r , unsigned max24bit )
{
	*q = (0x84210842ULL * max24bit)>>36;
	*r = max24bit + (*q) - ((*q)<<5);
	if(31<=*r) { *r-=31; *q+=1; } /* should we try to remove this ??? */
}

static inline uint32_t cvt_bin24_31x5( const uint8_t * gfv )
{
	unsigned r=gfv[4];
	r = (r<<5)-r+gfv[3];
	r = (r<<5)-r+gfv[2];
	r = (r<<5)-r+gfv[1];
	return (r<<5)-r+gfv[0];
}

static inline void cvt_31x5_bin24( uint8_t * gfv , uint32_t bv )
{
	uint32_t r;
	uint32_t bin24 = bv&0xffffff;
	div31( &bin24 , & r , bin24 ); gfv[0] = r;
	div31( &bin24 , & r , bin24 ); gfv[1] = r;
	div31( &bin24 , & r , bin24 ); gfv[2] = r;
	div31( &bin24 , & r , bin24 ); gfv[3] = r; gfv[4] = bin24;
}

static inline void cvt_bin96_31x20( uint32_t * bv , const uint8_t * gfv )
{
	uint32_t r;
	bv[0]=cvt_bin24_31x5(gfv);
	bv[1]=cvt_bin24_31x5(gfv+5);
	bv[2]=cvt_bin24_31x5(gfv+10);
	r    =cvt_bin24_31x5(gfv+15);
	bv[0] |= r<<24;
	bv[1] |= (r&0xff00)<<16;
	bv[2] |= (r&0xff0000)<<8;
}

static inline void cvt_31x20_bin96( uint8_t * gfv , const uint32_t * bv )
{
	uint32_t t = (bv[0]>>24)|((bv[1]>>16)&0xff00)|((bv[2]>>8)&0xff0000);
	cvt_31x5_bin24( gfv , bv[0] );
	cvt_31x5_bin24( gfv+5 , bv[1] );
	cvt_31x5_bin24( gfv+10 , bv[2] );
	cvt_31x5_bin24( gfv+15 , t );
}

static inline void pack_40b_31x8( uint8_t * pvec , const uint8_t * vec )
{
	uint32_t r = ((uint32_t)vec[0]) | (((uint32_t)vec[1])<<5) | (((uint32_t)vec[2])<<10) | (((uint32_t)vec[3])<<15) |
		(((uint32_t)vec[4])<<20) | (((uint32_t)vec[5])<<25) | (((uint32_t)vec[6])<<30);
	*((uint32_t *) (&pvec[0])) = r;
	pvec[4] = ((vec[7]<<3) | (vec[6]>>2));
}

static inline void unpack_31x8_40b( uint8_t * vec , const uint8_t *pvec )
{
	uint32_t ii;
	vec[0]=pvec[0]&0x1f;
	vec[1]=((*(const uint32_t*)pvec)>>5)&0x1f;
	ii = (*(const uint32_t*)(&pvec[1]))>>2;
	vec[2]=ii&0x1f; ii>>=5;
	vec[3]=ii&0x1f; ii>>=5;
	vec[4]=ii&0x1f; ii>>=5;
	vec[5]=ii&0x1f; ii>>=5;
	vec[6]=ii&0x1f; ii>>=5;
	vec[7]=ii&0x1f;
}

/* lnear algebra library */



static inline int vec_cmp40( const uint8_t * v1 , const uint8_t * v2 )
{
	int r = 0;
	unsigned i;
	for(i=0;i<40;i++) r |= (v1[i]^v2[i]);
	return r;
}


void vec_assign32( uint32_t * r , const uint8_t * inp , unsigned len );

void vec_assign( uint8_t * r , const uint8_t * inp , unsigned len );


void vec_fullreduce32_cvt( uint8_t * r , const uint32_t * inp , unsigned len );

void vec_fullreduce( uint8_t * r , unsigned len ); /* 0 - 62 ----> 0 - 30 */


void vec_setzero( uint8_t * vec , unsigned len );

void vec_negative( uint8_t * r , const uint8_t * frd_vec , unsigned len );

void vec_add( uint8_t *accu_r , const uint8_t *vec , unsigned len );


uint8_t vec_dot( const uint8_t * vec1 , const uint8_t * vec2 , unsigned len ); /* u */



void mat_mad32( uint32_t * accu_r , const uint8_t * mat , const uint8_t * vec , unsigned len );

void rowmat_mul32( uint8_t * r , const uint32_t * mat , const uint8_t * vec , unsigned len );


int solve_linear20( uint8_t * r , uint32_t * rowmat , const uint8_t * inp );



/* mq polynomial library */


#define NUM_QUAD_TERMS(n_var) ((n_var)*(n_var+1)/2)


typedef struct {
uint8_t l[64][40];
uint8_t q[NUM_QUAD_TERMS(64)][40];
} qpoly_64x40_t;


typedef struct {
uint8_t l[24][20];
uint8_t q[NUM_QUAD_TERMS(24)][20];
} qpoly_24x20_t;


typedef struct {
uint8_t l[44][20];
uint8_t q[NUM_QUAD_TERMS(44)][20];
} qpoly_44x20_t;



void eval_q64x40( uint8_t *r , const qpoly_64x40_t *poly , const uint8_t * inp );

void eval_q44x20( uint8_t *r , const qpoly_44x20_t *poly , const uint8_t * inp );

void eval_q24x20( uint8_t *r , const qpoly_24x20_t *poly , const uint8_t * inp );


void interpolate_64x40( qpoly_64x40_t * poly , void (*q_poly)(uint8_t *r,const void*,const uint8_t*), const void *key );



/* random functions */


void vec_rand( uint8_t * vec , unsigned len );

void mat_rand( uint8_t * mat , uint8_t * invmat , unsigned len ); /* len <= 64 , slow */




#if defined(__DEBUG__)

#include <stdio.h>

inline static void vec_dump( const char * s , const uint8_t * vec , unsigned len )
{
        unsigned i;
        printf("%s",s);
        for(i=0;i<len;i++) printf("%d, ",vec[i]);
        printf("\n");
}

inline static void vec_dump32( const char * s , const uint32_t * vec , unsigned len )
{
        unsigned i;
        printf("%s",s);
        for(i=0;i<len;i++) printf("%d, ",vec[i]);
        printf("\n");
}

inline static void mat_dump( const char * s , const uint8_t * vec , unsigned len )
{
        unsigned i;
        printf("%s\n",s);
        for(i=0;i<len;i++) vec_dump("",vec+len*i,len);
        printf("\n");
}

inline static void mat_dump32( const char * s , const uint32_t * vec , unsigned len )
{
        unsigned i;
        printf("%s\n",s);
        for(i=0;i<len;i++) vec_dump32("",vec+len*i,len);
        printf("\n");
}

#else

inline static void vec_dump( const char * s , const uint8_t * vec , unsigned len ) {}

inline static void vec_dump32( const char * s , const uint32_t * vec , unsigned len ) {}

inline static void mat_dump( const char * s , const uint8_t * vec , unsigned len ) {}

inline static void mat_dump32( const char * s , const uint32_t * vec , unsigned len ) {}

#endif /* defined(__DEBUG__) */

#ifdef __cplusplus
}
#endif

#endif /* _LINEAR_31_H_ */


