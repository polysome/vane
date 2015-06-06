#include <stdlib.h>

#define WORD	unsigned char
#define FIELD	256


#define AROW	8
#define ACOL	8  
#define BROW	8
#define BCOL	8
#define CROW	8
#define CCOL	8 
#define VARIABLE	64
#define EQUATION	128
#define CENTRAL_MAP_SIZE	2080
#define PUBLIC_KEY_SIZE		266240
#define SECRET_KEY_SIZE		294912

WORD add(WORD, WORD);
WORD sub(WORD, WORD); 
WORD mul(WORD, WORD);
WORD divs(WORD, WORD);
WORD inv(WORD);
int matrixinv(WORD *, int);
void matrixmul(WORD *, WORD *, int m, int n, int k, WORD *);
int matrixtranspose(WORD *,int,WORD *);
int tensorproduct(WORD *, WORD *, WORD *, int, int);
int echelonform(WORD *, int, int);
int keypair(WORD *, WORD *);
int encryption(WORD *, WORD *, WORD *, int);
int decryption(WORD *, WORD *, WORD *,  int);
