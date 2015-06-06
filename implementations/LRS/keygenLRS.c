//#include "keygen.h"
#include <stdio.h>
#include <stdlib.h>
#include "m4rie/m4rie.h"
#include "m4rie/gf2e.h"
#include "m4rie/mzed.h"
#include "m4rie/newton_john.h"

int counter, col_1, col_2;
FILE *fp;

int* gamma;
int** matrix_B1;
int** matrix_B2;

int** matrix_MS;
int** matrix_MS_inverse;
int** matrix_MS22_inverse;

int** matrix_MT;
int** matrix_MT_inverse;


int** matrix_A;
int** matrix_A_inverse_T;
int** matrix_A11_inverse_T;

int** matrix_Q11;
int** matrix_Q12;
int** matrix_Q21;
int** matrix_Q22;

int** matrix_F1;
int** matrix_F2;

int** matrix_MF;
int** matrix_MFK;
int** matrix_MPK;

int*** matrix_MQ;
int*** matrix_MFT;
int*** matrix_PK;

void free_gamma()
{
	free(gamma);
}

void free_matrix(int** matrix, size)
{
	int i;

	for(i = 0; i < size; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void free_3d_array(int*** matrix, int l, int w)
{
	int i, j;

	for(i = 0; i < l; i++)
	{
		for(j = 0; j < w; j++)
		{
			free(matrix[i][j]);
		}
		free(matrix[i]);
	}
	free(matrix);
}

void generate_gamma()
{
	gamma = malloc(m*sizeof(int));
	int i;

	for(i = 0; i < m; i++)
	{
		gamma[i] = (rand()%256);
	}
}
void generate_matrix_B1()
{
	matrix_B1 = malloc(m*(sizeof(int*)));
	int i;

	for(i = 0; i < m; i++)
	{
		matrix_B1[i] = malloc(D1*sizeof(int));
	}

	int j;

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < D1; j++)
		{
			if(j == 0)
			{
				matrix_B1[i][0] = 1;
			}
			else
			{
				matrix_B1[i][j] = mulTable(gamma[i], matrix_B1[i][j-1]);
			}
		}
	}	
}

void generate_matrix_MS()
{
	matrix_MS = malloc(m*sizeof(int*));
	int i, j;

	for(i = 0; i < m; i++)
	{
		matrix_MS[i] = malloc((m+1)*sizeof(int));
	}

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < (m+1); j++)
		{
			matrix_MS[i][j] = rand()%256;
		}
	}
}

void generate_matrix_MS_inverse()
{
	uint16_t int Irred[9]; /*Irreducible polynomials */
	Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
	Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
	Irred[8] = 1;

	gf2e *ff;
	int k;

	for(k = 2; k <= 8; k++)
	{
		ff = gf2e_init(Irred[k]);
	}

	const size_t mat_dim = m;
	int i,j;
	mzed_t *S = mzed_init(ff, ,mat_dim, mat_dim);

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < m; j++)
		{
			mzed_add_elem(S, i, j, matrix_MS[i][j]);
		}
	}
	mzed_t *S_INV = mzed_init(ff, mat_dim, mat_dim);
	mzed_invert_newton_john(S_INV, S);
	mzed_free(S);

	matrix_MS_inverse = malloc(m*sizeof(int*));

	for(i = 0; i < m; i++)
	{
		matrix_MS_inverse[i] = malloc((m+1)*sizeof(int));
	}

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < m; j++)
		{
			matrix_MS_inverse[i][j] = mzed_read_elem(S_INV, i, j);
		}
		matrix_MS_inverse[i][m] = matrix_MS[i][m];
	}
	mzed_free(S_INV);
	gf2e_free(ff);

}
void generate_matrix_MS22_inverse()
{
	uint16_t int Irred[9];
	Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
	Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
	Irred[8] = 1;

	gf2e *ff;
	int k;

	for(k = 2; k <= 8; k++)
	{
		ff = gf2e_init(Irred[k]);
	}

	const size_t mat_dim = o2;
	mzed_t *S22 = mzed_init(ff, mat_dim, mat_dim);
	int i, j;

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < o2; j++)
		{
			/* o1 + i >= o2 || o1 + j >= o2 .*/
			mzed_add_elem(S22, i, j, matrix_MS[o1+i][o1+j]);
		}
	}
	mzed_t *S22_INV = mzed_init(ff, mat_dim, mat_dim);
	mzed_invert_newton_john(S22_INV, S22);
	mzed_free(S22);

	matrix_MS22_inverse = malloc(o2*sizeof(int*));
	for(i = 0; i < o2; i++)
	{
		matrix_MS22_inverse[i] = malloc(o2*sizeof(int));
	}
	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < o2; j++)
		{
			matrix_MS22_inverse[i][j] = mzed_read_elem(S22_INV, i, j);
		}
	}
	mzed_free(S22_INV);
	gf2e_free(ff);
}

void generate_matrix_MT()
{
	matrix_MT = malloc((n+1)*sizeof(int*));

	int i, j;

	for(i = 0; i < n+1; i++)
	{
		matrix_MT[i] = malloc((n+1)*sizeof(int));
	}

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			matrix_MT[i][j] = rand()%256;
		}
	}

	for(i = 0; i < n; i++)
	{
		matrix_MT[n][i] = 0;
	}
	matrix_MT[n][n] = 1;
}

void generate_matrix_MT_inverse()
{
	uint16_t int Irred[9];
	Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
	Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
	Irred[8] = 1;

	gf2e *ff;
	int k;

	for(k = 2; k <= 8; k++)
	{
		ff = gf2e_init(Irred[k]);
	}

	const size_t mat_dim = n;
	mzed_t *T = mzed_init(ff, mat_dim, mat_dim);
	int i, j;

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			mzed_add_elem(T, i, j, matrix_MT[i][j]);
		}
	}

	mzed_t *T_INV = mzed_init(ff,mat_dim, mat_dim);
	mzed_invert_newton_john(T_INV, T);
	mzed_free(T);
	matrix_MT_inverse = malloc(n*sizeof(int*));

	for(i = 0; i < n; i++)
	{
		matrix_MT_inverse[i] = malloc((n+1)*int);
	}
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			matrix_MT_inverse[i][j] = mzed_read_elem(T_INV, i, j);
		}
		matrix_MT_inverse[i][n] = matrix_MT[i][n];
	}
	mzed_free(T_INV);
	gf2e_free(ff);
}

void generate_matrix_A()
{
	int r,s,i,j,p,q;
	matrix_A = malloc(D2*sizeof(int*));

	int k, l;

	for(k = 0; k < D2; k++)
	{
		matrix_A[k] = malloc(D2*sizeof(int));
	}

	i = j = p = 0;

	for(k = 0; k < D1; k++)
	{
		r = s = 0;

		for(l = 0; l < D1; l++)
		{
			if(i == j)
			{
				matrix_A[k][l] = mulTable(matrix_MT[r][i], matrix_MT[s][i]);
			}
			else
			{
				matrix_A[k][l] = addTable(mulTable(matrix_MT[r][i], matrix_MT[s][j]), mulTable(matrix_MT[s][i],matrix_MT[r][j]));
			}
			s++;

			if(s == v2)
			{
				r++;
				s=r;
			}
		}

		r = 0;
		s = v2;
		q = 0;

		for(l = D1; l < D2; l++)
		{
			if(i == j)
			{
				matrix_A[k][l]=mulTable(matrix_MT[r][i],matrix_MT[s][i]);
			}
			else
			{
				matrix_A[k][l]=addTable(mulTable(matrix_MT[r][i],matrix_MT[s][j]),mulTable(matrix_MT[s][i],matrix_MT[r][j]));
			}	

			if(r < v1)
			{
				s++;

				if(s==n)
				{
					r++;
					s = v2;
				}
			}

			if(r >= v1)
			{
				s = r + q;
				q++;

				if(s == n)
				{
					r++;
					s = r;
					q = 1;
				}
			}
		}
		j++;

		if(j == v2)
		{
			i++;
			j = i;
		}
	}

	i = 0;
	j = v2;

	for(k = D1; k < D2; k++)
	{
		r = s = 0;

		for(l = 0; l < D1; l++)
		{
			if(i == j)
			{
				matrix_A[k][l]=mulTable(matrix_MT[r][i],matrix_MT[s][i]);
			}
			else
			{
				matrix_A[k][l]=addTable(mulTable(matrix_MT[r][i],matrix_MT[s][j]),mulTable(matrix_MT[s][i],matrix_MT[r][j]));
			}
			s++;

			if(s == v2)
			{
				r++;
				s = r;
			}
		}

		r = 0;
		s = v2;
		q = 0;

		for(l = D1; l < D2; l++)
		{
			if( i == j)
			{
				matrix_A[k][l]=mulTable(matrix_MT[r][i],matrix_MT[s][i]);
			}
			else
			{
				matrix_A[k][l]=addTable(mulTable(matrix_MT[r][i],matrix_MT[s][j]),mulTable(matrix_MT[s][i],matrix_MT[r][j]));
			}

			if(r < v1)
			{
				s++;

				if(s == n)
				{
					r++;
					s = v2;
				}
			}

			if(r >= v1)
			{
				s = r + q;
				q++;

				if(s == n)
				{
					r++;
					s = r;
					q = 1;
				}
			}
		}

		if(i < v1)
		{
			j++;

			if(j == n)
			{
				i++;
				j = v2;
			}
		}

		if(i > = v1)
		{
			j = i+p;
			p++;

			if(j == n)
			{
				i++;
				j = i;
				p = 1;
			}
		}
	}
}

void generate_matrix_A_inverse()
{
	uint16_t int Irred[9];
	Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
	Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
	Irred[8] = 1;
	int k;
	gf2e *ff;

	for(k = 2; k <= 8; k++)
	{
		ff = gf2e_init(Irred[k]);
	}

	const size_t mat_dim = D2;
	mzed_t *A = mzed_init(ff, mat_dim, mat_dim);
	int i,j;

	for(i = 0; i < D2; i++)
	{
		for(j = 0; j < D2; j++)
		{
			mzed_add_elem(A, i, j, matrix_A[i][j]);
		}
	}
	mzed_t *A_INV = mzed_init(ff, mat_dim, mat_dim);
	mzed_invert_newton_john(A_INV, A);
	matrix_A_inverse_T = malloc(D2*sizeof(int*));
	mzed_free(A);

	for(i = 0; i < D2; i++)
	{
		matrix_A_inverse_T[i] = malloc(D2*sizeof(int));
	}

	for(i = 0; i < D2; i++)
	{
		for(j = 0; j < D2; j++)
		{
			matrix_A_inverse_T[j][i] = mzed_read_elem(A_INV, i, j);
		}
	}
	mzed_free(A_INV);
	gf2e_free(ff);
}

void generate_matrix_A11_inverse_T()
{
	uint16_t int Irred[9];
	Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
	Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
	Irred[8] = 1;
	int k;
	gf2e *ff;

	for(k = 2; k <= 8; k++)
	{
		ff = gf2e_init(Irred[k]);
	}

	const size_t mat_dim = D1;
	mzed_t *A11 = mzed_init(ff, mat_dim,mat_dim);
	int i, j;

	for(i = 0; i < D1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			mzed_add_elem(A11, i, j, matrix_A[i][j]);
		}
	}
	mzed_t *A11_INV = mzed_init(ff, mat_dim, mat_dim);
	mzed_invert_newton_john(A11_INV, A11);
	matrix_A_inverse_T = malloc(D1*sizeof(int*));
	mzed_free(A11);

	for(i = 0; i < D1; i++)
	{
		matrix_A_inverse_T[i] = malloc(D1*sizeof(int));
	}
	for(i = 0; i < D1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			matrix_A11_inverse_T[j][i] = mzed_read_elem(A11_INV, i, j);
		}
	}
	mzed_free(A11_INV);
	gf2e_free(ff);
}

void generate_matrix_Q11_Q21()
{
	matrix_Q11 = malloc(o1*sizeof(int*));
	int i, j, k;

	for(i = 0; i < o1; i++)
	{
		matrix_Q11[i] = malloc(D1*sizeof(D1));
	}

	matrix_Q21 = malloc(o2*sizeof(int*));

	for(i = 0; i < o2; i++)
	{
		matrix_Q21[i] = malloc(D1*sizeof(int));
	}

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			matrix_Q11[i][j] = 0;

			for(k = 0; k < m; k++)
			{
				matrix_Q11[i][j] = addTable(matrix_Q11[i][j], mulTable(matrix_MS_inverse[i][k], matrix_B1[k][j]));
			}
		}
	}

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < D2; j++)
		{
			matrix_Q21[i][j] = 0;

			for(k = 0; k < m; k++)
			{
				matrix_Q21[i][j] = addTable(matrix_Q21[i][j], mulTable(matrix_MS_inverse[o1+i][k], matrix_B1[k][j]));
			}
		}
	}	
}

void generate_matrix_Q12()
{
	matrix_Q12 = malloc(o1*sizeof(int*));

	int i, j, k;

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < D21; j++)
		{
			matrix_Q12[i][j] = 0;

			for(k = 0; k < D1; k++)
			{
				matrix_Q12[i][j] = addTable(matrix_Q12[i][j], mulTable(matrix_F1[i][k], matrix_A[D1+j][k]));
			}
		}
	}
}

void generate_matrix_Q22()
{
	matrix_Q22 = malloc(o2*sizeof(int*));

	int i, j;

	for(i = 0; i < o2; i++)
	{
		matrix_Q22[i] = malloc(D21*sizeof(int));
	}

	int** temp = malloc(o2*sizeof(int*));

	for(i = 0; i < o2; i++)
	{
		temp[i] = malloc(D21*sizeof(int));
	}

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < D21; j++)
		{
			temp[i][j] = 0;

			for(k = 0; k < o1; k++)
			{
				temp[i][j] = addTable(temp[i][j], mulTable(matrix_MS[o1+i][k], matrix_Q12[k][j]));
			}
		}
	}

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < D21; j++)
		{
			temp[i][j] = addTable(temp[i][j], matrix_B2[i][j]);
		}
	}

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < D21; j++)
		{
			matrix_Q22[i][j] = 0;

			for(k = 0; k < o2; k++)
			{
				matrix_Q22[i][j] = addTable(matrix_Q22[i][j], mulTable(matrix_MS22_inverse[i][k], temp[k][j]));
			}
		}
	}
}

void generate_matrix_F1()
{
	matrix_F1 = malloc(o1*sizeof(int*));

	int i, j , k;

	for(i = 0; i < o1; i++)
	{
		matrix_F1[i] = malloc(D1*sizeof(int));
	}

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			matrix_F1[i][j] = 0;

			for(k = 0; k < D1; k++)
			{
				matrix_F1[i][j] = addTable(matrix_F1[i][j], mulTable(matrix_Q11[i][k], matrix_A11_inverse_T[k][j]));
			}
		}
	}
}

void generate_matrix_F2()
{
	matrix_F2 = malloc(o2*sizeof(int*));

	int i, j, k;

	for(i = 0; i < o2; i++)
	{
		matrix_F2[i] = malloc(D2*sizeof(int));
	}

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < D2; j++)
		{
			matrix_F2[i][j] = 0;

			for(k = 0; k < D1; k++)
			{
				matrix_F2[i][j] = addTable(matrix_F2[i][j], mulTable(matrix_Q21[i][k], matrix_A_inverse_T[k][j]));
			}
			for(k = 0; k < D21; k++)
			{
				matrix_F2[i][j] = addTable(matrix_F2[i][j],mulTable( matrix_Q22[i][k],matrix_A_inverse_T[D1+k][j]));
			}
		}
	}
}

void generate_matrix_MF()
{
	col_1 = v1*o1 + (v1*(v1+1))/2 + v2 + 1;
	col_2 = v2*o2 + (v2*(v2+1))/2 + v2 + o2 + 1;

	matrix_MF = malloc(m*sizeof(int*));

	int i, j, k;

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			matrix_MF[i][j] = matrix_F1[i][j];
		}
		for(j = D1; j < col_1; j++)
		{
			matrix_MF[i][j] = rand()%256;
		}
	}

	for(i = 0; i < o2; i++)
	{
		int cf2 = 0;
		int cmf = 0;
		int cmf_2 = v2;

		for(j = 0; j < v1; j++)
		{
			for(k = j; k < v2; k++)
			{
				matrix_MF[o1+i][cmf] = matrix_F2[i][cf2];
				cf2++;
				cmf++;
			}
			cmf = cmf + o2;
		}

		for(j = 0; j < v1; j++)
		{
			for(k = 0; k < o2; k++)
			{
				matrix_MF[o1+i][cmf_2] = matrix_F2[i][cf2];
				cf2++;
				cmf_2++;
			}
			cmf_2 = cmf_2 + v2 - 1 - j;
		}

		for(j = D1+v1*o2; j<D2; j++)
		{
			matrix_MF[o1+i][j] = matrix_F2[i][j];
		}

		for(j = D2; j < col_2; j++)
		{
			matrix_MF[o1+i][j] = rand()%256;
		}
	}
}

void generate_privatekey()
{
	generate_gamma();
	/* seems like gamma is only needed for matrix B1 and B2. */
	generate_matrix_B1();
	generate_matrix_B2();
	/* free gamma right after generating B1 and B2. */
	free_gamma();

	generate_matrix_MS();
	generate_matrix_MS_inverse();
	generate_matrix_MS22_inverse();
	
	generate_matrix_MT();
	generate_matrix_MT_inverse();

	generate_matrix_A();
	generate_matrix_A_inverse_T();
	generate_matrix_A11_inverse_T();

	generate_matrix_Q11_Q21();
	generate_matrix_F1();
	generate_matrix_Q12();
	generate_matrix_Q22();
	generate_matrix_F2();
	
	generate_matrix_MF();

	/*Probably need to free matrix after generating, but should we output the key first.*/
}

void generate_matrix_MQ()
{
	matrix_MQ = malloc(m*sizeof(int**)); /*not sure*/
	int i,j,k, l;

	for(i = 0; i < m; i++)
	{
		matrix_MQ[i] = malloc((n+1)*sizeof(int*));

		for(j = 0; j < n+1; j++)
		{
			matrix_MQ[i][j] = malloc((n+1)*sizeof(int));
		}
	}

	for(i = 0; i < o1; i++)
	{
		counter = 0;

		for(k = 0; k < v1; k++)
		{
			for(l = 0; l < k; l++)
			{
				matrix_MQ[i][k][l] = 0;
			}

			for(l = k, l < v2; l++)
			{
				matrix_MQ[i][k][l] = matrix_MF[i][counter];
				counter++;
			}

			for(l = v2; l < n; l++)
			{
				matrix_MQ[i][k][l] = 0;
			}
		}

		for(k = v1; k < (n+1); k++)
		{
			for(l = 0; l < n; l++)
			{
				matrix_MQ[i][k][l] = 0;
			}
		}

		for(k = 0; k < v2; k++)
		{
			matrix_MQ[i][k][n] = matrix_MF[i][counter];
			counter++;
		}

		for(k = v2; k < n; k++)
		{
			matrix_MQ[i][k][n] = 0;
		}

		matrix_MQ[i][n][n] = matrix_MF[i][counter];
	}

	for(i = o1; i < m; i++)
	{
		counter = 0;

		for(k = 0; k < v2; k++)
		{
			for(l = 0; l < k; l++)
			{
				matrix_MQ[i][k][l] = 0;
			}
			for(l = k; l < n; l++)
			{
				matrix_MQ[i][k][l] = matrix_MF[i][counter];
				counter++;
			}
		}

		for(k = v2; k < (n+1); k++)
		{
			for(l = 0; l < n; l++)
			{
				matrix_MQ[i][k][l] = 0;
			}
		}

		for(k = 0; k < (n+1); k++)
		{
			matrix_MQ[i][k][n] = matrix_MF[i][counter];
			counter++;
		}
	}
}

void generate_matrix_MFT()
{
	int* temp = malloc((n+1)*sizeof(int));
	matrix_MFT = malloc(m*sizeof(int**));
	int i, j;

	for(i = 0; i < m; i++)
	{
		matrix_MFT[i] = malloc((n+1)*sizeof(int*));
		for(j = 0; j < n+1; j++)
		{
			matrix_MFT[i][j] = malloc((n+1)*sizeof(int));
		}
	}

	int l, k, q, r;

	for(i = 0; i < m; i++)
	{
		for(l = 0; l < n+1; l++)
		{
			for(k = 0; k < n+1; k++)
			{
				temp[k] = 0;

				for(p = 0; p < n+1; p++)
				{
					temp[m] = addTable(temp[m], mulTable(matrix_MT[p][l], matrix_MQ[i][p][k]));
				}
			}

			for(q = 0; q < n+1; q++)
			{
				matrix_MFT[i][l][q] = 0;

				for(r= 0; r < n+1; r++)
				{
					matrix_MFT[i][l][q]=addTable(matrix_MFT[i][l][q], mulTable(temp[r], matrix_MT[r][q]));
				}
			}
		}
	}
}

void generate_matrix_MFT_UT()
{
	int k, i, j;

	for(k = 0; k < m; k++)
	{
		for(i = 0; i < n+1; i++)
		{
			for(j = i+1; j < n+1; j++)
			{
				matrix_MFT[k][i][j]=addTable(matrix_MFT[k][i][j],matrix_MFT[k][j][i]);
				matrix_MFT[k][j][i]=0;
			}
		}
	}
}

void generate_matrix_MFK()
{
	matrix_MFK = malloc(m*sizeof(int*));
	int i, j, k, l;

	for(j = 0; j < m; j++)
	{
		matrix_MFK[j] = malloc(D*sizeof(int));
	}
	for(i = 0; i < m; i++)
	{
		counter = 0;

		for(k = 0; k < n; k++)
		{
			for(l = k; l < n; l++)
			{
				matrix_MFK[i][counter]=matrix_MFT[i][k][l];
				counter++;
			}
		}

		for(k = 0; k < n+1; k++)
		{
			matrix_MFK[i][counter]=matrix_MFT[i][k][n];
			counter++;
		}
	}
}

void generate_matrix_MPK()
{
	matrix_MPK = malloc(m*sizeof(int*));
	int i, j, k;

	for(j = 0; i < m; j++)
	{
		matrix_MPK[j] = malloc(D*sizeof(int));
	}

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < D; j++)
		{
			matrix_MPK[i][j] = 0;

			for(k = 0; k < m; k++)
			{
				matrix_MPK[i][j]=addTable(matrix_MPK[i][j],mulTable(matrix_MS[i][k],matrix_MFK[k][j]));
			}
		}
		matrix_MPK[i][D-1]=addTable(matrix_MPK[i][D-1],matrix_MS[i][m]);
	}
}

void generate_matrix_pk()
{
	matrix_PK = malloc(m*sizeof(int**));
	int i, j, k;

	for(i = 0; i < m; i++)
	{
		matrix_PK[i] = malloc((n+1)*sizeof(int*));

		for(j = 0; j < n+1; j++)
		{
			matrix_PK[i][j] = malloc((n+1)*sizeof(int));
		}
	}

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n; j++)
		{
			for(k = 0; k < j; k++)
			{
				matrix_PK[i][j][k] = 0;
			}

			for(k = j; k < n; k++)
			{
				matrix_PK[i][j][k]=matrix_MPK[i][j*n-(j*(j-1))/2+k-j];
			}
			matrix_PK[i][j][n]=matrix_MPK[i][(n*(n+1))/2 + j];
		}

		for(j = 0; j < n; j++)
		{
			matrix_p\[i][n][j] = 0;
		}
		matrix_PK[i][n][n]=matrix_MPK[i][((n+2)*(n+1))/2 -1];
	}
}

void generate_publickey()
{
	generate_matrix_MQ();
	generate_matrix_MFT();
	generate_matrix_MFT_UT();
	generate_matrix_MFK();
	generate_matrix_MPK();
	generate_matrix_pk();
}

void generate_key()
{
	generate_privatekey();
	generate_publickey();	
	/* Same issue here, when should free the matrix. */	
}

void print_info()
{
	fp = fopen("Key.txt", "w+");
	fprintf(fp, "%s", "\t\t\t********Key Generation********\n\n");
	fprintf(fp, "%s", "Instruction\n\n");
	fprintf(fp, "%s", "1. If matrix has more than 15 elements in a column than colums are divided into blocks of 15 elements\n   and new elements start from next line.\n\n");
	fprintf(fp, "%s", "2. New row always start one line below the previous row.\n\n");
	fprintf(fp, "%s", "3. Private Key in Matrix(MS,MT and MF) and Public Key in both Matrix(MPK) and Polynomial form are printed.\n\n\n");
	fflush(fp);
	fclose(fp);
}

void print_matrix_B1()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix B2 is:\n\n");
	int i, j;

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < D21; j++)
		{
			fprintf(fp, "\t%d", matrix_B2[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_MS()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix MS including vector cs  is:\n\n");
	int i, j;

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < m+1; j++)
		{
			fprintf(fp, "\t%d", matrix_MS[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_MS_inverse()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix MS inverse is:\n\n");
	int i, j;

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < m; j++)
		{
			fprintf(fp, "\t%d", matrix_MS_inverse[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_MS22_inverse()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix S22 inverse is:\n\n");
	int i, j;

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < o2; j++)
		{
			fprintf(fp, "\t%d", matrix_MS22_inverse[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_MT()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix MT including vector ct is:\n\n");
	int i, j;

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n+1; j++)
		{
			fprintf(fp, "\t%d", matrix_MT[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_MT_inverse()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix MT inverse is:\n\n");
	int i, j;

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			fprintf(fp, "\t%d", matrix_MT_inverse[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_A()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix A is:\n\n");
	int i, j;

	for(i = 0; i < D2; i++)
	{
		for(j = 0; j < D2; j++)
		{
			fprintf(fp, "\t%d", matrix_A[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_A_inverse_T()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix A inverse transposed is:\n\n");
	int i, j;

	for(i = 0; i < D2; i++)
	{
		for(j = 0; j < D2; j++)
		{
			fprintf(fp, "\t%d", matrix_A_inverse_T[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_A11_inverse_T()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix A11 inverse transposed is:\n\n");
	int i, j;

	for(i = 0; i < D1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			fprintf(fp, "\t%d", matrix_A11_inverse_T[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_Q11_Q21()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix Q11 is:\n\n");
	int i, j;

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			fprintf(fp, "\t%d", matrix_Q11[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fprintf(fp, "%s", "\n\nMatrix Q21 is:\n\n");

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			fprintf(fp, "\t%d", matrix_Q21[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_F1()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix F1 is:\n\n");
	int i, j;

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < D1; j++)
		{
			fprintf(fp, "\t%d", matrix_F1[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_Q12()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix Q12 is:\n\n");
	int i, j;

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < D21; j++)
		{
			fprintf(fp, "\t%d", matrix_Q12[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_Q22()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix Q22 is:\n\n");
	int i, j;

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < D21; j++)
		{
			fprintf(fp, "\t%d", matrix_Q22[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_F2()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix F2 is:\n\n");
	int i, j;

	for(i = 0; i < o2; i++)
	{
		for(j = 0; j < D2; j++)
		{
			fprintf(fp, "\t%d", matrix_F2[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_MF()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix MF is:\n\n");
	int i, j;

	for(i = 0; i < o1; i++)
	{
		for(j = 0; j < (v1*o1 + (v1*(v1+1))/2 + v2 + 1); j++)
		{
			fprintf(fp, "\t%d", matrix_MF[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);

	for(i = o1; i < m; i++)
	{
		for(j = 0; j < (v2*o2 + (v2*(v2+1))/2 + v2 + o2 + 1); j++)
		{
			fprintf(fp, "\t%d", matrix_MF[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
		}
		fflush(fp);
	}
	fflush(fp);
	fclose(fp);
}

void print_matrix_MPK()
{
	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nMatrix MPK is:\n\n");
	int i, j;

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < ((n+1)*(n+2))/2; j++)
		{
			fprintf(fp, "\t%d", matrix_MPK[i][j]);
			if((j+1)%15 == 0)
			{
				fprintf(fp, "%s\n", "");
			}
			fflush(fp);
		}
	}
	fflush(fp);
	fclose(fp);
}

void print_key()
{
	print_info();
	print_matrix_B1();
	print_matrix_B2();

	print_matrix_MS();
	print_matrix_MT();
	print_matrix_MF();
	print_matrix_MPK();

	fp = fopen("Key.txt", "a+");
	fprintf(fp, "%s", "\n\nPublic Key is written in polynomial form\n\n");
	int i;

	for(i = 0; i < m; i++)
	{
		fprintf(fp, "%s%d%s", "P(", i+1, "):\t");

		int k = 1;
		int l = 1;
		int j;

		for(j = 0; j < ((n)*(n+1))/2; j++)
		{
			fprintf(fp, "%d%s%d%s%d\n", );
			if(matrix_MPK[i][j] <= 9 && k < 10 && l <10)
			{
				fprintf(fp, "%s\t", "");
				fflush(fp);
			}
			fprintf(fp, "%s\t", "+");

			if(l%n == 0)
			{
				k++;
				l = k;
			}
			else
			{
				l++;
			}
			if((j+1)%5 == 0)
			{
				fprintf(fp, "%s\t", "");
				fflush(fp);
			}
		}
		fprintf(fp, "%s\t", "");

		for(j = 0; j < n; j++)
		{
			fprintf(fp, "%d%s%d\t\t+\t", matrix_MPK[i][(n*(n+1))/2 + j], "*m", j+1);
			if((j+1)%5 == 0)
			{
				fprintf(fp, "%s\t", "");
				fflush(fp);
			}
		}
		fprintf(fp, "%s\t", "");
		fprintf(fp, "%d\n", matrix_MPK[i][((n+1)*(n+2))/2 -1]);
		fflush(fp);
	}
	fclose(fp);
}

void store_privatekey()
{
	fp = fopen("PrivateKeycyclic.txt", "w+");
	fprintf(fp, "%s","o1\n");
	fprintf(fp, "%d\n", o1);

	fprintf(fp, "%s","o2\n");
	fprintf(fp, "%d\n", o2);

	fprintf(fp, "%s","v1\n");
	fprintf(fp, "%d\n", v1);

	fprintf(fp, "%s","v2\n");
	fprintf(fp, "%d\n", v2);

	fprintf(fp, "%s","m\n");
	fprintf(fp, "%d\n", m);

	fprintf(fp, "%s","n\n");
	fprintf(fp, "%d\n", n);

	fprintf(fp, "%s","MF\n");
	int k, j;
	for(k = 0; k < o1; k++)
	{
		for(j = 0; j < col_1; j++)
		{
			fprintf(fp, "%d%s", matrix_MF[k][j], " \n");
		}
		//fprintf(fp, "%s\n","");
		fflush(fp);
	}
	for(k = o1; k < m; k++)
	{
		for(j = 0; j < col_2; j++)
		{
			fprintf(fp, "%d%s", matrix_MF[k][j], " \n");
		}
		//fprintf(fp, "%s\n","");
		fflush(fp);
	}
	fprintf(fp, "%s\n", "S");

	for(k = 0; k < m; k++)
	{
		for(j = 0; j < m+1; j++)
		{
			fprintf(fp, "%d%s", matrix_MS_inverse[k][j], " \n");
		}
		//fprintf(fp, "%s\n","");
		fflush(fp);
	}
	fprintf(fp, "%s\n", "T");

	for(k = 0; k < n; k++)
	{
		for(j = 0; j < n+1; j++)
		{
			fprintf(fp, "%d%s", matrix_MS_inverse[k][j], " \n");
		}
		//fprintf(fp, "%s\n","");
		fflush(fp);
	}
	fclose(fp);
}

void store_publickey()
{
	fp = fopen("PublicKeycyclic.txt", "w+");

	fprintf(fp, "%s","o1\n");
	fprintf(fp, "%d\n", o1);

	fprintf(fp, "%s","o2\n");
	fprintf(fp, "%d\n", o2);

	fprintf(fp, "%s","v1\n");
	fprintf(fp, "%d\n", v1);

	fprintf(fp, "%s","v2\n");
	fprintf(fp, "%d\n", v2);

	fprintf(fp, "%s","m\n");
	fprintf(fp, "%d\n", m);

	fprintf(fp, "%s","n\n");
	fprintf(fp, "%d\n", n);

	fprintf(fp, "%s\n", "pk");
	int i,j,k;

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n+1; j++)
		{
			fprintf(fp, "%d ", matrix_PK[i][j][k]);
		}
		fprintf(fp, "%s\n","");
		fflush(fp);
	}
	fclose(fp);
}


