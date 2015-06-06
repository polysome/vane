#ifndef PUBLICKEYGEN_H
#define PUBLICKEYGEN_H

#include "basic.h"
#include <m4rie/finite_field_givaro.h>
#include <m4rie/m4rie.h>
#include <m4rie/gf2e_matrix.h>


extern int** matrix_MS;
extern int** matrix_MT;
extern int** matrix_MF;
extern int** matrix_MS_inverse;
extern int** matrix_MT_inverse;

extern int*** matrix_PK;

void generate_matrix_MS();
void generate_matrix_MT();//Generates Matrix MT of size (n+1)*(n+1) by using coefficients and constants of Linear Map T 
void generate_matrix_MF();//It generates matrix MF
void generate_privatekey();//It generates polynomial equations

void generate_matrix_MQ();//Generates Matrix MQ of size (n+1)*(n+1) by using coefficients and constants of Quadratic Polynomial
void generate_matrix_MFT();//Generates matrix MP=Transpose(MT)*MQ*MT
void generate_matrix_MFT_UT();//Generates matrix MP in triangular form
void generate_matrix_MFK();
void generate_matrix_MPK();//Generates matrix MPK
void generate_matrix_pk();//Generates Matrix PK of size (n+1)*(n+1) by using Matrix MPK(Public Key)
void generate_publickey();

void generate_matrix_MS_inverse();
void generate_matrix_MT_inverse();
void generate_Linear_map_inverse();

void generate_key();

void print_info();//It prints the basic Instruction in Text file.
void print_matrix_MS();
void print_matrix_MT();//It prints Matrix MT in Project_A.txt
void print_matrix_MF();//It prints Matrix MF in Project_A.txt
void print_matrix_MPK();//It prints Matrix MPK in Project_A.txt
void print_key();//Prints Public Key in Polynomial form

void store_privatekey();
void store_publickey();




#endif
