#ifndef BASIC_H
#define BASIC_H



//#include <ctime> 
//#include "keygen.h"


 //extern int v1;
 //extern int v2;
 //extern int o1;
 //extern int o2;
 //extern int m;
 //extern int n;//Oilvariables+VinegarVariables
 //extern int D;  //number of quadratic terms in a public polynomial
 //extern int D1,D2,D21;

 int addTable(int,int);//It add two elemnts using Galois Arithmetic
 int mulTable(int,int);//It multiplies two elements using Galois Arithmetic for GF(2^8) and Primitive Polynomial x^8 + x^4 + x^3 + x + 1 
 int* divisionlookuparray();//It stores Multiplicative Inverse of each elemnt
 int divideTable(int,int,int []);//It gives the division result of two elements by using lookup array
 int singlerandomvaluegenerator();//It generates a single random value
 void randvectorgenerator(int,int*);//It generates vector with random values
 void randvalformatrix(int,int,int**);//It generates matrix with random elements using row,columns and matrix (in which elements are stored) as input
 unsigned long long rdtsc(void);//It counts the number of CPU cycles
 double* ave_max_min(double [],int);//Returns the Maximum,Average and Minimum of Time 
 unsigned long long* ave_max_min_cpu(unsigned long long [],int);//Returns the Maximum,Average and Minimum of Cpu Cycles



#endif