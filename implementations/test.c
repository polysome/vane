#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "abc.c"


int 
main(int argc, char *argv[])
{
	int plaintextlen = VARIABLE;
	int i;
	WORD sk[SECRET_KEY_SIZE];
	WORD pk[PUBLIC_KEY_SIZE];
	WORD plaintext[VARIABLE];
	WORD ciphertext[EQUATION];
	WORD decrypttext[EQUATION];
	srand((int)time(0));
    keypair(sk, pk);
	printf( "\n plaintext: \n");

	for (i = 0; i < VARIABLE; i++) {
		plaintext[i] = rand() % FIELD ;
		printf( "%d ", plaintext[i]) ;
	}
	printf("\n");
	encryption(pk, plaintext, ciphertext, plaintextlen);
	printf( "\n ciphertext: ");
	for (i = 0; i < EQUATION; i++) {
		printf( "%d ", ciphertext[i]) ;
	}
	printf("\n");
	decryption(sk, decrypttext, ciphertext,  EQUATION);
	printf( "\n decrypttext: \n");
	for (i = 0; i < VARIABLE; i++) {
		printf( "%d ", decrypttext[i]) ;
	}
	printf("\n");
	return (1);
} 