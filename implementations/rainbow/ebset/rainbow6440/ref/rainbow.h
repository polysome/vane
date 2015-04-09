
#ifndef _RAINBOW_H_
#define _RAINBOW_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "linear31.h"


typedef qpoly_64x40_t pubkey_t;


typedef struct {
uint8_t s[64][64];
uint8_t sc[64];
qpoly_24x20_t vv1st;
uint8_t ov1st_rowmat[20][20][24];
uint8_t ol1st_rowmat[20][20];
qpoly_44x20_t vv2nd;
uint8_t ov2nd_rowmat[20][20][44];
uint8_t ol2nd_rowmat[20][20];
uint8_t t[40][40];
uint8_t tc[40];
} seckey_t;


/* functions on GF(31) */

#if defined(__DEBUG__)
int genkey( uint8_t * pubkey , uint8_t * seckey );
int sign( uint8_t * s , const seckey_t * key , const uint8_t * m );
int verify( const uint8_t * md , const pubkey_t * key , const uint8_t * s );
#endif

/* public interface */

#define SECKEY_SIZE_BYTE 60960
#define PUBKEY_SIZE_BYTE 53600
#define DIGEST_SIZE_BYTE 24
#define SIGNATURE_SIZE_BYTE 40


int genkey_pack( uint8_t * pubkey , uint8_t *seckey );

int sign_bin( uint8_t * s320b , const seckey_t * key , const uint8_t * md192b );

int verify_bin( const uint8_t * md192b , const uint8_t * key , const uint8_t * s320b );


#ifdef __cplusplus
}
#endif

#endif /* _RAIDBOW_H_ */



