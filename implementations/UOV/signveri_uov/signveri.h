#ifndef SIGNVERI_H
#define SIGNVERI_H


extern int result;
//extern int* hash_verify;

#include "basic.h"




void generate_matrix_pk(int);//Generates Matrix PK of size (n+1)*(n+1) by using Matrix MPK(Public Key)
void signverification();//Generates Hash values by Hash_verify=(signature,1)*PK*transpose(signature,1) and compare with original hashvalue for authentication process

void print_signatureverification(int);//It prints hashvalues generated by verification process and result of authentication.





#endif
