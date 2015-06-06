#ifndef KEYGEN_H
#define KEYGEN_H

//#include "../m4rie/m4rie.h"

//#include "basic.h"

extern int* Gamma;
extern int** matrix_B1;
extern int** matrix_B2;

extern int** matrix_MS;
extern int** matrix_MS_inverse;
extern int** matrix_MS22_inverse;

extern int** matrix_MT;
extern int** matrix_MT_inverse;


extern int** matrix_A;
extern int** matrix_A_inverse_T;
extern int** matrix_A11_inverse_T;

extern int** matrix_Q11;
extern int** matrix_Q12;
extern int** matrix_Q21;
extern int** matrix_Q22;

extern int** matrix_F1;
extern int** matrix_F2;

extern int** matrix_MF;
extern int** matrix_MFK;
extern int** matrix_MPK;

extern int*** matrix_MQ;
extern int*** matrix_MFT;
extern int*** matrix_PK;

extern int D;
extern int D1;
extern int D2;
extern int D21;
extern int m, n;
extern int o1, o2, v1, v2;




void free_gamma();
void free_matrix(int** matrix, int size);
void free_3d_array(int*** matrix, int l, int w);
void generate_gama();
void generate_matrix_B1();
void generate_matrix_B2();
void generate_matrix_MS();
void generate_matrix_MS_inverse();
void generate_matrix_MS22_inverse();
void generate_matrix_MT();
void generate_matrix_MT_inverse();
void generate_matrix_A();
void generate_matrix_A_inverse_T();
void generate_matrix_A11_inverse_T();
void generate_matrix_Q11_Q21();
void generate_matrix_Q12();
void generate_matrix_Q22();
void generate_matrix_F1();
void generate_matrix_F2();
void generate_matrix_MF();
void generate_privatekey();
void generate_matrix_MQ();
void generate_matrix_MFT();
void generate_matrix_MFT_UT();
void generate_matrix_MFK();
void generate_matrix_MPK();
void generate_matrix_pk();
void generate_publickey();
void generate_key();
void print_info();
void print_matrix_B1();
void print_matrix_B2();
void print_matrix_MS();
void print_matrix_MS_inverse();
void print_matrix_MS22_inverse();
void print_matrix_MT();
void print_matrix_MT_inverse();
void print_matrix_A();
void print_matrix_A_inverse_T();
void print_matrix_A11_inverse_T();
void print_matrix_Q11_Q21();
void print_matrix_F1();
void print_matrix_Q12();
void print_matrix_Q22();
void print_matrix_F2();
void print_matrix_MF();
void print_matrix_MPK();
void print_key();
void store_privatekey();
void store_publickey();
void free_all_privatekey();
void free_all_publickey();


#endif