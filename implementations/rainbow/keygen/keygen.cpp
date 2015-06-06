#include "keygen.h"

using namespace std;
using namespace M4RIE;

ofstream myfile_pk;


int counter;
int col_1,col_2;

int** matrix_MS;
int** matrix_MT;
int** matrix_MS_inverse;
int** matrix_MT_inverse;
int** matrix_MF;

int*** matrix_MQ;
int*** matrix_MFT;
int** matrix_MFK;
int** matrix_MPK;
int*** matrix_PK;

void generate_matrix_MS()
{	
	

	matrix_MS = new int *[m] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < (m) ; i++ )
		matrix_MS[i] = new int[m+1];//memory allocated for  elements of each column.
	
	for(int i=0;i<m;i++)
		for(int j=0;j<(m+1);j++)
			matrix_MS[i][j]=rand()%256;
	
	
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



void generate_matrix_MF()
{	
	matrix_MF = new int *[m] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o1 ; i++ )
	matrix_MF[i] = new int[v1*o1 + (v1*(v1+1))/2 + v2 + 1];//memory allocated for  elements of each column.
	
	for( int i = o1 ; i < m ; i++ )
	matrix_MF[i] = new int[v2*o2 + (v2*(v2+1))/2 + v2 + o2 + 1];//memory allocated for  elements of each column.
	
	col_1=v1*o1 + (v1*(v1+1))/2 + v2 + 1,col_2=v2*o2 + (v2*(v2+1))/2 + v2 + o2 + 1;
	
	for(int i=0;i<o1;i++)
		for(int j=0;j<col_1;j++)
			matrix_MF[i][j]=rand()%256;

	for(int i=o1;i<m;i++)
		for(int j=0;j<col_2;j++)
			matrix_MF[i][j]=rand()%256;
					
}
void generate_privatekey()
{	
	generate_matrix_MS();
	generate_matrix_MT();
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
void generate_Linear_map_inverse()
{
	generate_matrix_MS_inverse();
	generate_matrix_MT_inverse();		
}
void generate_key()
{
	generate_privatekey();
	generate_publickey();	
	generate_Linear_map_inverse();
	
}


void print_info()
{	
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_conventional_output/Key.txt");
	myfile_pk<<"\t\t\t********Key Generation********\n\n";
	myfile_pk<<"Instruction\n\n";
	myfile_pk<<"1. If matrix has more than 15 elements in a column than colums are divided into blocks of 15 elements\n   and new elements start from next line.\n\n";
	myfile_pk<<"2. New row always start one line below the previous row.\n\n";
	myfile_pk<<"3. Private Key in Matrix(MS,MT and MF) and Public Key in both Matrix(MPK) and Polynomial form are printed.\n\n\n";
	myfile_pk.close();

}
void print_matrix_MS()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_conventional_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MS is:\n\n";

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
void print_matrix_MT()
{		
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_conventional_output/Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MT is:\n\n";

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

void print_matrix_MF()
{
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_conventional_output/Key.txt",ios::app);
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
	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_conventional_output/Key.txt",ios::app);
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
	print_matrix_MS();
	print_matrix_MT();
	print_matrix_MF();
	print_matrix_MPK();


	myfile_pk.open("/home/ishtiaq/rainbow/rainbow_conventional_output/Key.txt",ios::app);
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
	myfile_pk.open("/home/ishtiaq/rainbow/data/PrivateKey.txt");
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
	myfile_pk.open("/home/ishtiaq/rainbow/data/PublicKey.txt");
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
