#ifndef SIGN_H
#define SIGN_H 

#include "basic.h"

extern int* hashvalue;
extern int* signature;
 
void hashvaluegenerator();//It generates hashvalues
void randomvalues4vinegarVariables();//It generates random values for vinegar variables 
void linearmap_of_S(); 
void Quadtraticmap();
void linearmap_of_T();

int* guasselimination(int**,int,int);//It returns solution of matrix by taking number of rows,number of columns and matrix as input
 
void generatesignature();//It generates signature of message 
 
void print_sign_generation();
void store_hash_sign();


#endif
