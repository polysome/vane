#include "keygen.h"

using namespace std;
using namespace M4RIE;

ofstream myfile_pk;


int counter,col_1,col_2;

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

void generate_matrix_B1()
{	
	matrix_B1 = new int *[m] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < (m) ; i++ )
		matrix_B1[i] = new int[D1];//memory allocated for  elements of each column.
	
	for(int i=0;i<D1;i++)
		matrix_B1[0][i]=rand()%256;

	 for(int i=1;i<m;i++)
	 {
		 for(int j=0;j<D1-i;j++)   //Enteringelemnts to upper triangular matrix
		 {
			 matrix_B1[i][j+i]=matrix_B1[0][j];

		 }
		 for(int k=D1;k>D1-i;k--)//Entering elemnts to lower triangular matrix
		 {
			 matrix_B1[i][k-D1+i-1]=matrix_B1[0][k-1];
		 }
	 }

}
void generate_matrix_B2()
{
	matrix_B2 = new int *[o2] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < (o2) ; i++ )
		matrix_B2[i] = new int[D21];//memory allocated for  elements of each column.
	
	for(int i=0;i<D21;i++)
		matrix_B2[0][i]=rand()%256;

	 for(int i=1;i<o2;i++)
	 {
		 for(int j=0;j<D21-i;j++)   //Enteringelemnts to upper triangular matrix
		 {
			 matrix_B2[i][j+i]=matrix_B2[0][j];

		 }
		 for(int k=D21;k>D21-i;k--)//Entering elemnts to lower triangular matrix
		 {
			 matrix_B2[i][k-D21+i-1]=matrix_B2[0][k-1];
		 }
	 }
}
void generate_matrix_MS()
{	
	

	matrix_MS = new int *[m] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < (m) ; i++ )
		matrix_MS[i] = new int[m+1];//memory allocated for  elements of each column.
	
	for(int i=0;i<m;i++)
		for(int j=0;j<(m+1);j++)
			matrix_MS[i][j]=rand()%256;
	
	
}
void generate_matrix_MS_inverse()
{
	std::vector< GFqDom<long>::Residu_t > Irred(9);
		Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
		Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
		Irred[8] = 1;

		FiniteField *F = (FiniteField*)(new GFqDom<long>(2,8,Irred));
		gf2e *ff = gf2e_init_givgfq(F);
	
	const size_t mat_dim=m;
		
	mzed_t *S = mzed_init(ff,mat_dim,mat_dim);
		
	for(int i=0;i<m;i++)
		for(int j=0;j<m;j++)
			mzed_add_elem(S,i,j,matrix_MS[i][j]);
	
	mzed_t *S_INV = mzed_init(ff,mat_dim,mat_dim);

	mzed_invert_travolta(S_INV,S); 
		
	matrix_MS_inverse = new int *[m] ;//memory allocated for elements of rows.
		for( int i = 0 ; i < m ; i++ )
			matrix_MS_inverse[i] = new int[m+1];//memory allocated for  elements of each column.
		
		for(int i=0;i<m;i++)
		{
			for(int j=0;j<m;j++)
				matrix_MS_inverse[i][j]=mzed_read_elem(S_INV,i,j);

			matrix_MS_inverse[i][m]=matrix_MS[i][m];
		}
			

}
void generate_matrix_MS22_inverse()
{
	std::vector< GFqDom<long>::Residu_t > Irred(9);
		Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
		Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
		Irred[8] = 1;

		FiniteField *F = (FiniteField*)(new GFqDom<long>(2,8,Irred));
		gf2e *ff = gf2e_init_givgfq(F);
	
	const size_t mat_dim=o2;
		
	mzed_t *S22 = mzed_init(ff,mat_dim,mat_dim);
		
	for(int i=0;i<o2;i++)
		for(int j=0;j<o2;j++)
			mzed_add_elem(S22,i,j,matrix_MS[o1+i][o1+j]);
	
	mzed_t *S22_INV = mzed_init(ff,mat_dim,mat_dim);

	mzed_invert_travolta(S22_INV,S22); 
		
	matrix_MS22_inverse = new int *[o2] ;//memory allocated for elements of rows.
		for( int i = 0 ; i < o2 ; i++ )
			matrix_MS22_inverse[i] = new int[o2];//memory allocated for  elements of each column.
		
		for(int i=0;i<o2;i++)
			for(int j=0;j<o2;j++)
				matrix_MS22_inverse[i][j]=mzed_read_elem(S22_INV,i,j);


}


void generate_matrix_MT()
{	
	

	matrix_MT = new int *[n+1] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < (n+1) ; i++ )
		matrix_MT[i] = new int[n+1];//memory allocated for  elements of each column.
	
	for(int i=0;i<n;i++)
		for(int j=0;j<(n+1);j++)
			matrix_MT[i][j]=rand()%256;
	
	
	for(int i=0;i<n;i++)
		matrix_MT[n][i]=0;

	matrix_MT[n][n]=1;

}
void generate_matrix_MT_inverse()
{
	std::vector< GFqDom<long>::Residu_t > Irred(9);
		Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
		Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
		Irred[8] = 1;

		FiniteField *F = (FiniteField*)(new GFqDom<long>(2,8,Irred));
		gf2e *ff = gf2e_init_givgfq(F);
	
	const size_t mat_dim=n;
		
	mzed_t *T = mzed_init(ff,mat_dim,mat_dim);
		
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			mzed_add_elem(T,i,j,matrix_MT[i][j]);
	
	mzed_t *T_INV = mzed_init(ff,mat_dim,mat_dim);

	mzed_invert_travolta(T_INV,T); 
		
	matrix_MT_inverse = new int *[n] ;//memory allocated for elements of rows.
		for( int i = 0 ; i < n ; i++ )
			matrix_MT_inverse[i] = new int[n+1];//memory allocated for  elements of each column.
		
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
				matrix_MT_inverse[i][j]=mzed_read_elem(T_INV,i,j);

			matrix_MT_inverse[i][n]=matrix_MT[i][n];
		}
}
void generate_matrix_A()
{	
	int r,s,i,j,p,q;
	matrix_A = new int *[D2] ;//memory allocated for elements of rows.
	
	for( int k = 0 ; k < D2 ; k++ )

	matrix_A[k] = new int[D2];//memory allocated for  elements of each column.
	
	i=0;
	j=0;
	p=0;
	
	for(int k=0;k<D1;k++)
	{	
		r=0;
		s=0;
			
		for(int l=0;l<D1;l++)
		{	
			//cout<<i+1<<j+1<<r+1<<s+1<<"\t";
			if(i==j)
				matrix_A[k][l]=mulTable(matrix_MT[r][i],matrix_MT[s][i]);
			else
				matrix_A[k][l]=addTable(mulTable(matrix_MT[r][i],matrix_MT[s][j]),mulTable(matrix_MT[s][i],matrix_MT[r][j]));	
			s++;
			
			if(s==v2)
			{
				r++;
				s=r;
			}
			
			
			
		}

		r=0;
		s=v2;
		q=0;
		for(int l=D1;l<D2;l++)
		{	
			//cout<<i+1<<j+1<<r+1<<s+1<<"\t";
			if(i==j)
				matrix_A[k][l]=mulTable(matrix_MT[r][i],matrix_MT[s][i]);
			else
				matrix_A[k][l]=addTable(mulTable(matrix_MT[r][i],matrix_MT[s][j]),mulTable(matrix_MT[s][i],matrix_MT[r][j]));	
			if(r<v1)
			{
				s++;
		
				if(s==n)
				{
				  r++;
				  s=v2;
				}
			}
				
			if(r>=v1)
			{	
				s=r+q;
				q++;				
				if(s==n)
				{
					r++;
					s=r;
					q=1;
				}
						
			}						

			
		}
		j++;
			
			if(j==v2)
			{
				i++;
				j=i;
			}
		//cout<<"\n";


	}
	
	i=0;
	j=v2;
	
	for(int k=D1;k<D2;k++)
	{	
		r=0;
		s=0;
			
		for(int l=0;l<D1;l++)
		{	
			//cout<<i+1<<j+1<<r+1<<s+1<<"\t";
			if(i==j)
				matrix_A[k][l]=mulTable(matrix_MT[r][i],matrix_MT[s][i]);
			else
				matrix_A[k][l]=addTable(mulTable(matrix_MT[r][i],matrix_MT[s][j]),mulTable(matrix_MT[s][i],matrix_MT[r][j]));	
			s++;
			
			if(s==v2)
			{
				r++;
				s=r;
			}
			
			
			
		}

		r=0;
		s=v2;
		q=0;
		for(int l=D1;l<D2;l++)
		{	
			//cout<<i+1<<j+1<<r+1<<s+1<<"\t";
			if(i==j)
				matrix_A[k][l]=mulTable(matrix_MT[r][i],matrix_MT[s][i]);
			else
				matrix_A[k][l]=addTable(mulTable(matrix_MT[r][i],matrix_MT[s][j]),mulTable(matrix_MT[s][i],matrix_MT[r][j]));	
			if(r<v1)
			{
				s++;
		
				if(s==n)
				{
				  r++;
				  s=v2;
				}
			}
				
			if(r>=v1)
			{	
				s=r+q;
				q++;				
				if(s==n)
				{
					r++;
					s=r;
					q=1;
				}
						
			}						

			
		}
			
		if(i<v1)
		{
			j++;
		
			if(j==n)
			{
			  i++;
			  j=v2;
			}
		}
				
		if(i>=v1)
		{	
			j=i+p;
			p++;				
			if(j==n)
			{
				i++;
				j=i;
				p=1;
			}
						
		}	
	//cout<<"\n";
	}
	
}
void generate_matrix_A_inverse_T()
{
	std::vector< GFqDom<long>::Residu_t > Irred(9);
		Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
		Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
		Irred[8] = 1;

		FiniteField *F = (FiniteField*)(new GFqDom<long>(2,8,Irred));
		gf2e *ff = gf2e_init_givgfq(F);
	
	const size_t mat_dim=D2;
		
	mzed_t *A = mzed_init(ff,mat_dim,mat_dim);
		
	for(int i=0;i<D2;i++)
		for(int j=0;j<D2;j++)
			mzed_add_elem(A,i,j,matrix_A[i][j]);
	
	mzed_t *A_INV = mzed_init(ff,mat_dim,mat_dim);

	mzed_invert_travolta(A_INV,A); 
		
	matrix_A_inverse_T = new int *[D2] ;//memory allocated for elements of rows.
		for( int i = 0 ; i < D2 ; i++ )
			matrix_A_inverse_T[i] = new int[D2];//memory allocated for  elements of each column.
		
		for(int i=0;i<D2;i++)
			for(int j=0;j<D2;j++)
				matrix_A_inverse_T[j][i]=mzed_read_elem(A_INV,i,j);

}
void generate_matrix_A11_inverse_T()
{
	std::vector< GFqDom<long>::Residu_t > Irred(9);
		Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
		Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
		Irred[8] = 1;

		FiniteField *F = (FiniteField*)(new GFqDom<long>(2,8,Irred));
		gf2e *ff = gf2e_init_givgfq(F);
	
	const size_t mat_dim=D1;
		
	mzed_t *A11 = mzed_init(ff,mat_dim,mat_dim);
		
	for(int i=0;i<D1;i++)
		for(int j=0;j<D1;j++)
			mzed_add_elem(A11,i,j,matrix_A[i][j]);
	
	mzed_t *A11_INV = mzed_init(ff,mat_dim,mat_dim);

	mzed_invert_travolta(A11_INV,A11); 
		
	matrix_A11_inverse_T = new int *[D1] ;//memory allocated for elements of rows.
		for( int i = 0 ; i < D1 ; i++ )
			matrix_A11_inverse_T[i] = new int[D1];//memory allocated for  elements of each column.
		
		for(int i=0;i<D1;i++)
			for(int j=0;j<D1;j++)
				matrix_A11_inverse_T[j][i]=mzed_read_elem(A11_INV,i,j);

}
void generate_matrix_Q11_Q21()
{	
	matrix_Q11 = new int *[o1] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o1 ; i++ )
		matrix_Q11[i] = new int[D1];//memory allocated for  elements of each column.
	
	matrix_Q21 = new int *[o2] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o2 ; i++ )
		matrix_Q21[i] = new int[D1];//memory allocated for  elements of each column.

	for(int i=0;i<o1;i++)
		for(int j=0;j<D1;j++)
		{
			matrix_Q11[i][j]=0;
			
			for(int k=0;k<m;k++)
				matrix_Q11[i][j]=addTable(matrix_Q11[i][j],mulTable( matrix_MS_inverse[i][k],matrix_B1[k][j]));
		}
	for(int i=0;i<o2;i++)
		for(int j=0;j<D1;j++)
		{
			matrix_Q21[i][j]=0;
			
			for(int k=0;k<m;k++)
				matrix_Q21[i][j]=addTable(matrix_Q21[i][j],mulTable( matrix_MS_inverse[o1+i][k],matrix_B1[k][j]));
		}

}
void generate_matrix_Q12()
{
	matrix_Q12 = new int *[o1] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o1 ; i++ )
		matrix_Q12[i] = new int[D21];//memory allocated for  elements of each column.

	for(int i=0;i<o1;i++)
		for(int j=0;j<D21;j++)
		{
			matrix_Q12[i][j]=0;
			
			for(int k=0;k<D1;k++)
				matrix_Q12[i][j]=addTable(matrix_Q12[i][j],mulTable( matrix_F1[i][k],matrix_A[D1+j][k]));
		}

	
}
void generate_matrix_Q22()
{
	matrix_Q22 = new int *[o2] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o2 ; i++ )
		matrix_Q22[i] = new int[D21];//memory allocated for  elements of each column.

	int** temp= new int *[o2] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o2 ; i++ )
		temp[i] = new int[D21];//memory allocated for  elements of each column.

	for(int i=0;i<o2;i++)
		for(int j=0;j<D21;j++)
		{
			temp[i][j]=0;
			
			for(int k=0;k<o1;k++)
				temp[i][j]=addTable(temp[i][j],mulTable( matrix_MS[o1+i][k],matrix_Q12[k][j]));
		}
	for(int i=0;i<o2;i++)
		for(int j=0;j<D21;j++)
			temp[i][j]=addTable(temp[i][j],matrix_B2[i][j]);

	
	for(int i=0;i<o2;i++)
		for(int j=0;j<D21;j++)
		{
			matrix_Q22[i][j]=0;
			
			for(int k=0;k<o2;k++)
				matrix_Q22[i][j]=addTable(matrix_Q22[i][j],mulTable( matrix_MS22_inverse[i][k],temp[k][j]));
		}

	
	
}
void generate_matrix_F1()
{
	matrix_F1= new int *[o1] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o1 ; i++ )
		matrix_F1[i] = new int[D1];//memory allocated for  elements of each column.

	for(int i=0;i<o1;i++)
		for(int j=0;j<D1;j++)
		{
			matrix_F1[i][j]=0;
			
			for(int k=0;k<D1;k++)
				matrix_F1[i][j]=addTable(matrix_F1[i][j],mulTable( matrix_Q11[i][k],matrix_A11_inverse_T[k][j]));
		}
	
	
	
}
void generate_matrix_F2()
{	
	matrix_F2= new int *[o2] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o2 ; i++ )
		matrix_F2[i] = new int[D2];//memory allocated for  elements of each column.

	for(int i=0;i<o2;i++)
		for(int j=0;j<D2;j++)
		{
			matrix_F2[i][j]=0;
			
			for(int k=0;k<D1;k++)
				matrix_F2[i][j]=addTable(matrix_F2[i][j],mulTable( matrix_Q21[i][k],matrix_A_inverse_T[k][j]));
	
			for(int k=0;k<D21;k++)
				matrix_F2[i][j]=addTable(matrix_F2[i][j],mulTable( matrix_Q22[i][k],matrix_A_inverse_T[D1+k][j]));
		}

}

void generate_matrix_MF()
{	
	col_1=v1*o1 + (v1*(v1+1))/2 + v2 + 1,col_2=v2*o2 + (v2*(v2+1))/2 + v2 + o2 + 1;	
	
	matrix_MF = new int *[m] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o1 ; i++ )
	matrix_MF[i] = new int[col_1];//memory allocated for  elements of each column.
	
	for( int i = o1 ; i < m ; i++ )
	matrix_MF[i] = new int[col_2];//memory allocated for  elements of each column.
	
	
	
	for(int i=0;i<o1;i++)
	{
		for(int j=0;j<D1;j++)			
			matrix_MF[i][j]=matrix_F1[i][j];
			
		for(int j=D1;j<col_1;j++)
			matrix_MF[i][j]=rand()%256;		
	}
	
	for(int i=0;i<o2;i++)
	{	
		int cf2=0,cmf=0,cmf_2=v2;
		for(int j=0;j<v1;j++)
		{
			for(int k=j;k<v2;k++)
			{
				matrix_MF[o1+i][cmf]=matrix_F2[i][cf2];
				cf2++;
				cmf++;
			}
			cmf=cmf+o2;
		}
		for(int j=0;j<v1;j++)
		{
			for(int k=0;k<o2;k++)
			{
				matrix_MF[o1+i][cmf_2]=matrix_F2[i][cf2];
				cf2++;
				cmf_2++;
			}
			cmf_2=cmf_2+v2-1-j;
		}
		
	        for(int j=D1+v1*o2;j<D2;j++)
			matrix_MF[o1+i][j]=matrix_F2[i][j];
			
		
		for(int j=D2;j<col_2;j++)
			matrix_MF[o1+i][j]=rand()%256;		
	}
	
					
}
void generate_privatekey()
{	
	generate_matrix_B1();
	generate_matrix_B2();

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
		
}
void generate_matrix_MQ()
{
	
	matrix_MQ = new int **[m];
  	for (int i = 0; i < m; i++) {
    		matrix_MQ[i] = new int*[n+1];
		for (int j = 0; j < (n+1); j++)
      			matrix_MQ[i][j] = new int[n+1];
		
	}
	
	for(int i=0;i<o1;i++)
	{

		counter=0;
	
		for(int k=0;k<v1;k++)
		{	
			for(int l=0;l<k;l++)
				matrix_MQ[i][k][l]=0;
	
			for(int l=k;l<v2;l++)
			{
				matrix_MQ[i][k][l]=matrix_MF[i][counter];
				counter=counter+1;
			}
			
			for(int l=v2;l<n;l++)
				matrix_MQ[i][k][l]=0;
		}
		
		for(int k=v1;k<(n+1);k++)
			for(int l=0;l<n;l++)
				matrix_MQ[i][k][l]=0;
			

		for(int k=0;k<v2;k++)
		{
			matrix_MQ[i][k][n]=matrix_MF[i][counter];
			counter=counter+1;
		}

		for(int k=v2;k<n;k++)
			matrix_MQ[i][k][n]=0;
	
		matrix_MQ[i][n][n]=matrix_MF[i][counter];
	}
	
	for(int i=o1;i<m;i++)
	{

		counter=0;
	
		for(int k=0;k<v2;k++)
		{	
			for(int l=0;l<k;l++)
				matrix_MQ[i][k][l]=0;
	
			for(int l=k;l<n;l++)
			{
				matrix_MQ[i][k][l]=matrix_MF[i][counter];
				counter=counter+1;
			}
			
		}
		
		for(int k=v2;k<(n+1);k++)
			for(int l=0;l<n;l++)
				matrix_MQ[i][k][l]=0;
			

		for(int k=0;k<(n+1);k++)
		{
			matrix_MQ[i][k][n]=matrix_MF[i][counter];
			counter=counter+1;
		}

	
	}
	
}

void generate_matrix_MFT()
{
	int* temp= new int [n+1];
	
	matrix_MFT = new int **[m];
  	for (int i = 0; i < m; i++) {
    		matrix_MFT[i] = new int*[n+1];
		for (int j = 0; j < (n+1); j++)
      			matrix_MFT[i][j] = new int[n+1];
		
	}

	for(int i=0;i<m;i++)
	{	
		for(int l=0;l<(n+1);l++)
		{	
			for(int m=0;m<(n+1);m++)
			{	
				temp[m]=0;

				for(int p=0;p<(n+1);p++)
					temp[m]=addTable(temp[m],mulTable(matrix_MT[p][l],matrix_MQ[i][p][m]));

			}

			for(int q=0;q<(n+1);q++)
			{
				matrix_MFT[i][l][q]=0;

				for(int r=0;r<(n+1);r++)
					matrix_MFT[i][l][q]=addTable(matrix_MFT[i][l][q],mulTable(temp[r],matrix_MT[r][q]));
			}

		}
	}

}
void generate_matrix_MFT_UT()
{	
	for(int k=0;k<m;k++)
	{
	for(int i=0;i<(n+1);i++)
		for(int j=i+1;j<(n+1);j++)
		{
			matrix_MFT[k][i][j]=addTable(matrix_MFT[k][i][j],matrix_MFT[k][j][i]);
			matrix_MFT[k][j][i]=0;
		}
	}
}
void generate_matrix_MFK()
{
	matrix_MFK = new int *[m] ;//memory allocated for elements of rows.
	
	for( int j = 0 ; j < m ; j++ )
		matrix_MFK[j] = new int[D];//memory allocated for  elements of each column.
	
	for(int i=0;i<m;i++)
	{
		counter=0;

		for(int k=0;k<n;k++)
			for(int l=k;l<n;l++)
			{
				matrix_MFK[i][counter]=matrix_MFT[i][k][l];
				counter=counter+1;
			}

		for(int k=0;k<(n+1);k++)
		{
			matrix_MFK[i][counter]=matrix_MFT[i][k][n];
			counter=counter+1;
		}
	}
	

}
void generate_matrix_MPK()
{	
	matrix_MPK = new int *[m] ;//memory allocated for elements of rows.
	
	for( int j = 0 ; j < m ; j++ )
		matrix_MPK[j] = new int[D];//memory allocated for  elements of each column.
	
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<D;j++)
		{	
			matrix_MPK[i][j]=0;
			for(int k=0;k<m;k++)
			matrix_MPK[i][j]=addTable(matrix_MPK[i][j],mulTable(matrix_MS[i][k],matrix_MFK[k][j]));
		}
		matrix_MPK[i][D-1]=addTable(matrix_MPK[i][D-1],matrix_MS[i][m]);
	}
	
	
}
void generate_matrix_pk()
{
	matrix_PK = new int **[m];
  	for (int i = 0; i < m; i++) {
    		matrix_PK[i] = new int*[n+1];
		for (int j = 0; j < (n+1); j++)
      			matrix_PK[i][j] = new int[n+1];
	}
	
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			for(int k=0;k<j;k++)
				matrix_PK[i][j][k]=0;

			for(int k=j;k<n;k++)
				matrix_PK[i][j][k]=matrix_MPK[i][j*n-(j*(j-1))/2+k-j];

			matrix_PK[i][j][n]=matrix_MPK[i][(n*(n+1))/2 + j];

		}
	
		for(int j=0;j<n;j++)
			matrix_PK[i][n][j]=0;

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
		
}


void print_info()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt");
	myfile_pk<<"\t\t\t********Key Generation********\n\n";
	myfile_pk<<"Instruction\n\n";
	myfile_pk<<"1. If matrix has more than 15 elements in a column than colums are divided into blocks of 15 elements\n   and new elements start from next line.\n\n";
	myfile_pk<<"2. New row always start one line below the previous row.\n\n";
	myfile_pk<<"3. Private Key in Matrix(MS,MT and MF) and Public Key in both Matrix(MPK) and Polynomial form are printed.\n\n\n";
	myfile_pk.close();

}
void print_matrix_B1()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix B1 is:\n\n";

	for(int i=0;i<(m);i++)
	{
		for(int j=0;j<(D1);j++)
		{
			myfile_pk<<"\t"<<matrix_B1[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_B2()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix B2 is:\n\n";

	for(int i=0;i<(o2);i++)
	{
		for(int j=0;j<(D21);j++)
		{
			myfile_pk<<"\t"<<matrix_B2[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_MS()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MS including vector cs  is:\n\n";

	for(int i=0;i<(m);i++)
	{
		for(int j=0;j<(m+1);j++)
		{
			myfile_pk<<"\t"<<matrix_MS[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_MS_inverse()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MS inverse is:\n\n";

	for(int i=0;i<(m);i++)
	{
		for(int j=0;j<(m);j++)
		{
			myfile_pk<<"\t"<<matrix_MS_inverse[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}

void print_matrix_MS22_inverse()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix S22 inverse is:\n\n";

	for(int i=0;i<o2;i++)
	{
		for(int j=0;j<o2;j++)
		{
			myfile_pk<<"\t"<<matrix_MS22_inverse[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}

void print_matrix_MT()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MT including vector ct is:\n\n";

	for(int i=0;i<(n);i++)
	{
		for(int j=0;j<(n+1);j++)
		{
			myfile_pk<<"\t"<<matrix_MT[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_MT_inverse()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MT inverse is:\n\n";

	for(int i=0;i<(n);i++)
	{
		for(int j=0;j<(n);j++)
		{
			myfile_pk<<"\t"<<matrix_MT_inverse[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}

void print_matrix_A()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix A is:\n\n";

	for(int i=0;i<D2;i++)
	{
		for(int j=0;j<D2;j++)
		{
			myfile_pk<<"\t"<<matrix_A[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_A_inverse_T()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix A inverse transposed is:\n\n";

	for(int i=0;i<D2;i++)
	{
		for(int j=0;j<D2;j++)
		{
			myfile_pk<<"\t"<<matrix_A_inverse_T[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_A11_inverse_T()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix A11 inverse transposed is:\n\n";

	for(int i=0;i<D1;i++)
	{
		for(int j=0;j<D1;j++)
		{
			myfile_pk<<"\t"<<matrix_A11_inverse_T[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_Q11_Q21()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix Q11 is:\n\n";

	for(int i=0;i<o1;i++)
	{
		for(int j=0;j<D1;j++)
		{
			myfile_pk<<"\t"<<matrix_Q11[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;
	}
	myfile_pk<<"\n\nMatrix Q21 is:\n\n";

	for(int i=0;i<o1;i++)
	{
		for(int j=0;j<D1;j++)
		{
			myfile_pk<<"\t"<<matrix_Q21[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_F1()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix F1 is:\n\n";

	for(int i=0;i<o1;i++)
	{
		for(int j=0;j<D1;j++)
		{
			myfile_pk<<"\t"<<matrix_F1[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}
void print_matrix_Q12()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix Q12 is:\n\n";

	for(int i=0;i<o1;i++)
	{
		for(int j=0;j<D21;j++)
		{
			myfile_pk<<"\t"<<matrix_Q12[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}	
void print_matrix_Q22()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix Q22 is:\n\n";

	for(int i=0;i<o2;i++)
	{
		for(int j=0;j<D21;j++)
		{
			myfile_pk<<"\t"<<matrix_Q22[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}	
void print_matrix_F2()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix F2 is:\n\n";

	for(int i=0;i<o2;i++)
	{
		for(int j=0;j<D2;j++)
		{
			myfile_pk<<"\t"<<matrix_F2[i][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();

}	


void print_matrix_MF()
{
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MF  is:\n\n";

	for(int k=0;k<o1;k++)
	{
		for(int j=0;j<(v1*o1 + (v1*(v1+1))/2 + v2 + 1);j++)
		{
			myfile_pk<<"\t"<<matrix_MF[k][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}

		myfile_pk<<endl<<endl;

	}
	for(int k=o1;k<m;k++)
	{
		for(int j=0;j<(v2*o2 + (v2*(v2+1))/2 + v2 + o2 + 1);j++)
		{
			myfile_pk<<"\t"<<matrix_MF[k][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	
	myfile_pk.close();

}
void print_matrix_MPK()
{
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MPK is:\n\n";

	for(int k=0;k<m;k++)
	{
		for(int j=0;j<((n+1)*(n+2))/2;j++)
		{
			myfile_pk<<"\t"<<matrix_MPK[k][j];
			if((j+1)%15==0)
				myfile_pk<<"\n";
		}
		myfile_pk<<endl<<endl;

	}
	myfile_pk.close();
}
void print_key()
{	
	print_info();
	print_matrix_B1();
	print_matrix_B2();

	print_matrix_MS();
	//print_matrix_MS_inverse();
	//print_matrix_MS22_inverse();
	
	print_matrix_MT();
	//print_matrix_MT_inverse();

	//print_matrix_A();
	//print_matrix_A_inverse_T();
	//print_matrix_A11_inverse_T();

	//print_matrix_Q11_Q21();
	//print_matrix_F1();
	//print_matrix_Q12();
	//print_matrix_Q22();
	//print_matrix_F2();

	print_matrix_MF();
	print_matrix_MPK();


	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Key.txt",ios::app);
	myfile_pk<<"\n\nPublic Key is written in polynomail form\n\n";
	
	for(int i=0;i<m;i++)
	{	
		myfile_pk<<endl;
		myfile_pk<<"P("<<i+1<<"):\t";

		int k=1;
		int l=1;

		for(int j=0;j<((n)*(n+1))/2;j++)
		{
			myfile_pk<<matrix_MPK[i][j]<<"*m"<<k<<"*m"<<l<<"\t";
			if (matrix_MPK[i][j]<=9 && k<10 && l<10)
			myfile_pk<<"\t";
					
			myfile_pk<<"+\t";

			if(l%n==0)
			{	
				k++;
				l=k;
				
			}
			else
				l++;
					
			if((j+1)%5==0)
				myfile_pk<<endl<<"\t";
					
		}
		myfile_pk<<endl<<"\t";

		for (int j=0;j<n;j++)
		{
			myfile_pk<<matrix_MPK[i][(n*(n+1))/2 + j]<<"*m"<<j+1<<"\t\t+\t";
			if((j+1)%5==0)
				myfile_pk<<endl<<"\t";
		}
		

		myfile_pk<<endl<<"\t";
	
		myfile_pk<<matrix_MPK[i][((n+1)*(n+2))/2 -1]<<endl;
	}
	myfile_pk.close();

}
void store_privatekey()
{
	myfile_pk.open("/home/ishtiaq/rainbow/data/PrivateKeycyclic.txt");
	myfile_pk<<"o1\n";
	myfile_pk<<o1<<"\n";

	myfile_pk<<"o2\n";
	myfile_pk<<o2<<"\n";

	myfile_pk<<"v1\n";
	myfile_pk<<v1<<"\n";

	myfile_pk<<"v2\n";
	myfile_pk<<v2<<"\n";

	myfile_pk<<"m\n";
	myfile_pk<<m<<"\n";

	myfile_pk<<"n\n";
	myfile_pk<<n<<"\n";
	
	
	myfile_pk<<"MF\n";
	for(int k=0;k<o1;k++)
	{
		for(int j=0;j<col_1;j++)
			myfile_pk<<matrix_MF[k][j]<<" ";
	        myfile_pk<<endl;
	}
	for(int k=o1;k<m;k++)
	{
		for(int j=0;j<col_2;j++)
			myfile_pk<<matrix_MF[k][j]<<" ";
	        myfile_pk<<endl;
	}
	
	myfile_pk<<"S\n";
	for(int k=0;k<m;k++)
	{
		for(int j=0;j<(m+1);j++)
			myfile_pk<<matrix_MS_inverse[k][j]<<" ";
	        myfile_pk<<endl;
	}

	myfile_pk<<"T\n";
	for(int k=0;k<n;k++)
	{
		for(int j=0;j<(n+1);j++)
			myfile_pk<<matrix_MT_inverse[k][j]<<" ";

	        myfile_pk<<endl;
	}
	
	myfile_pk.close();
			
		
}
void store_publickey()
{
	myfile_pk.open("/home/ishtiaq/rainbow/data/PublicKeycyclic.txt");
	
	myfile_pk<<"o1\n";
	myfile_pk<<o1<<"\n";

	myfile_pk<<"o2\n";
	myfile_pk<<o2<<"\n";

	myfile_pk<<"v1\n";
	myfile_pk<<v1<<"\n";

	myfile_pk<<"v2\n";
	myfile_pk<<v2<<"\n";

	myfile_pk<<"m\n";
	myfile_pk<<m<<"\n";
	
	myfile_pk<<"n\n";
	myfile_pk<<n<<"\n";
	
	myfile_pk<<"pk\n";
	
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<(n+1);j++)
		{
			for(int k=0;k<(n+1);k++)
			myfile_pk<<matrix_PK[i][j][k]<<" ";
		myfile_pk<<endl;
		}
	}
	myfile_pk.close();
	
	

}
