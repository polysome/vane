#ifndef _CONFIG_H_
#define _CONFIG_H_

#define QUARTZ103



#ifdef QUARTZ103
#define CORE_SIZE 103
#define MAX_DEG 129
#define VINEGAR 4
#define MINUS 3
#define REPEAT 4

#define PUBKEY_BYTES 75514
#define SECKEY_BYTES 3774

#endif


#define N (CORE_SIZE+VINEGAR)
#define M (CORE_SIZE-MINUS)

#define SECMSG_BYTES (((M+(MINUS+VINEGAR)*REPEAT)+7)/8)


#define S_WIDTH N
#define T_WIDTH CORE_SIZE

#define PUBKEY_NUM_TERMS ((N)*(N+1)/2)

#define SIZE_BYTE_M ((M+7)>>3)
#define SIZE_BYTE_N ((N+7)>>3)
#define SIZE_BYTE_CORE ((CORE_SIZE+7)>>3)




#endif
