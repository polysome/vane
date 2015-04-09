#ifndef PUBLICKEYGEN_H
#define PUBLICKEYGEN_H

#include "basic.h"
#include <m4rie/finite_field_givaro.h>
#include <m4rie/m4rie.h>
#include <m4rie/gf2e_matrix.h>

extern int** coeffa;
extern int** coeffb;
extern int** coeffc;
extern int** coeffd;
extern int* coeffe;
extern int** coefft;
extern int* coeffv;

extern int** matrix_MPK;

void coeffagenerator();//It generates coefficents a matrix (secret coefficients)
void coeffbgenerator();//It generates coefficents b matrix (secret coefficients)
void coeffcgenerator();//It generates coefficents c matrix (secret coefficients)
void coeffdgenerator();//It generates coefficents d matrix (secret coefficients)
void coeffegenerator();//It generates coefficents e vector (secret coefficients)
void generatepolyeqs();//It generates polynomial equations
void coefftgenerator();//It generates coefficents t matrix (secret coefficients)
void coeffvgenerator();//It generates coefficents v vector (secret coefficients)

void generate_matrix_MT();//Generates Matrix MT of size (n+1)*(n+1) by using coefficients and constants of Linear Map T 
void generate_L_inverse();//Generates inverse of Linear Map
void generate_matrix_MQ(int i);//Generates Matrix MQ of size (n+1)*(n+1) by using coefficients and constants of Quadratic Polynomial
void generate_matrix_MF(int i);//It generates matrix MF
void generate_matrix_MP();//Generates matrix MP=Transpose(MT)*MQ*MT
void generate_matrix_MP_UT();//Generates matrix MP in triangular form
void generate_matrix_MPK();//Generates matrix MPK

void print_info();//It prints the basic Instruction in Text file.
void print_matrix_MT();//It prints Matrix MT in Project_A.txt
void print_matrix_MF();//It prints Matrix MF in Project_A.txt
void print_matrix_MPK();//It prints Matrix MPK in Project_A.txt
void print_publickey();//Prints Public Key in Polynomial form
void print_privatekey();//Prints Private Key in Polynomial form
void store_privatekey();
void store_publickey();





#endif
