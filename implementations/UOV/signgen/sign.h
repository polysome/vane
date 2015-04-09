#ifndef SIGN_H
#define SIGN_H 

#include "basic.h"


 
 extern int* hashvalue;
 extern int* message;

 int* randomvalues4vinegarVariables();//It generates random values for vinegar variables
 void hashvaluegenerator();//It generates hashvalues
 int* solvepolynomialequations();//It solves polynomial equations and returns the solution
 int* linearmapT(int*);//It returns signature of message by constructing linear map,it takes oil and vinegar variables as input
 int* guasselimination(int**,int,int);//It returns solution of matrix by taking number of rows,number of columns and matrix as input
 
 void generatesignature();//It generates signature of message 
 
 void print_sign_generation();
 void store_hash_sign();


#endif
