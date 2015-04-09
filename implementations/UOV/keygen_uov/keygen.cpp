#include "keygen.h"

using namespace std;
using namespace M4RIE;

ofstream myfile_pk;


int** coeffa;
int** coeffb;
int** coeffc;
int** coeffd;
int* coeffe;
int** coefft;
int** coefft_inv;
int* coeffv;

int** matrix_MT;
int** matrix_MQ;
int** matrix_MF;
int** matrix_MP;
int** matrix_MPK;


// coefficients of the central map ============================================================
void coeffagenerator()
{
	coeffa = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < oilvariables ; i++ )
		coeffa[i] = new int[oilvariables*vinegarvariables];//memory allocated for  elements of each column.
	randvalformatrix(oilvariables,oilvariables*vinegarvariables,coeffa);
}

void coeffbgenerator()
{
	coeffb = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < oilvariables ; i++ )
		coeffb[i] = new int[(vinegarvariables*(vinegarvariables+1))/2];//memory allocated for  elements of each column.
	randvalformatrix(oilvariables,(vinegarvariables*(vinegarvariables+1))/2,coeffb);
}

void coeffcgenerator()
{
	coeffc = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < oilvariables ; i++ )
		coeffc[i] = new int[vinegarvariables];//memory allocated for  elements of each column.
	randvalformatrix(oilvariables,vinegarvariables,coeffc);
}

void coeffdgenerator()
{
	coeffd = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < oilvariables ; i++ )
		coeffd[i] = new int[oilvariables];//memory allocated for  elements of each column.
	randvalformatrix(oilvariables,oilvariables,coeffd);
}


void coeffegenerator()
{
	coeffe = new int [oilvariables] ;
	randvectorgenerator(oilvariables,coeffe);

}
void generatepolyeqs()
{
	coeffagenerator();
	coeffbgenerator();
	coeffcgenerator();
	coeffdgenerator();
	coeffegenerator();
	
}
void coefftgenerator()
{
	coefft = new int *[oilvariables+vinegarvariables] ;//memory allocated for elements of rows.
		
	for( int i = 0 ; i < oilvariables+vinegarvariables ; i++ )//memory allocated for  elements of each column.
		coefft[i] = new int[oilvariables+vinegarvariables];
	
	randvalformatrix(oilvariables+vinegarvariables,oilvariables+vinegarvariables,coefft);
}

void coeffvgenerator()
{
	coeffv = new int [oilvariables+vinegarvariables];//memory allocated for elements of vector.
	randvectorgenerator(oilvariables+vinegarvariables,coeffv);
}

// affine map MT =============================================================
void generate_matrix_MT()
{	
	coefftgenerator();
	coeffvgenerator();

	matrix_MT = new int *[n+1] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < (n+1) ; i++ )
		matrix_MT[i] = new int[n+1];//memory allocated for  elements of each column.
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
			matrix_MT[i][j]=coefft[i][j];

		matrix_MT[i][n]=coeffv[i];
	}
	
	for(int i=0;i<n;i++)
		matrix_MT[n][i]=0;

	matrix_MT[n][n]=1;

}

// inversion ================================================================================
void generate_L_inverse()
{	
	std::vector< GFqDom<long>::Residu_t > Irred(9);
		Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
		Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
		Irred[8] = 1;

		FiniteField *F = (FiniteField*)(new GFqDom<long>(2,8,Irred));
		gf2e *ff = gf2e_init_givgfq(F);
	
	
	
	

		const size_t m=n;
		
		
		mzed_t *T = mzed_init(ff,m,m);
		

		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
				mzed_add_elem(T,i,j,coefft[i][j]);
		
                
		mzed_t *T_INV = mzed_init(ff,m,m);
  		mzed_invert_travolta(T_INV,T); 
		
		coefft_inv = new int *[n] ;//memory allocated for elements of rows.
		for( int i = 0 ; i < n ; i++ )
			coefft_inv[i] = new int[n];//memory allocated for  elements of each column.
		
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
				coefft_inv[i][j]=mzed_read_elem(T_INV,i,j);


}

// matrix MQ representing central polynomials (one upper triangular matrix for each polynomial)
void generate_matrix_MQ(int i)
{
	matrix_MQ = new int *[n+1] ;//memory allocated for elements of rows.
	
	for( int l = 0 ; l < (n+1) ; l++ )
		matrix_MQ[l] = new int[n+1];//memory allocated for  elements of each column.
	
	
	
	for(int j=0;j<vinegarvariables;j++)
	{	
		for(int l=0;l<j;l++)
			matrix_MQ[j][l]=0;

		for(int l=j;l<vinegarvariables;l++)
			matrix_MQ[j][l]=coeffb[i][vinegarvariables*j+l-j-((j*(j-1))/2)];

		for(int k=0;k<oilvariables;k++)
			matrix_MQ[j][k+vinegarvariables]=coeffa[i][k*vinegarvariables+j];

		matrix_MQ[j][n]=coeffc[i][j];
		
	}

	for(int j=0;j<oilvariables;j++)
	{
		for(int k=0;k<n;k++)
			matrix_MQ[j+vinegarvariables][k]=0;

		matrix_MQ[j+vinegarvariables][n]=coeffd[i][j];
	}

	for(int l=0;l<n;l++)
		matrix_MQ[n][l]=0;

	matrix_MQ[n][n]=coeffe[i];


}
// Maccauley matrix of the central map =========================
void generate_matrix_MF(int i)
{	
	matrix_MF[i] = new int[(vinegarvariables*(n+oilvariables+1)+2*n+2)/2];//memory allocated for  elements of each column.

	int counter=0;
	
	for(int k=0;k<vinegarvariables;k++)
		for(int l=k;l<n;l++)
		{
			matrix_MF[i][counter]=matrix_MQ[k][l];
			counter=counter+1;
		}
	for(int k=0;k<(n+1);k++)
	{
		matrix_MF[i][counter]=matrix_MQ[k][n];
		counter=counter+1;
	}
}

// matrix MP representing the public key (one matrix for each polynomial)
void generate_matrix_MP()
{
	matrix_MP = new int *[n+1] ;//memory allocated for elements of rows.
	
	for( int j = 0 ; j < (n+1) ; j++ )
		matrix_MP[j] = new int[n+1];//memory allocated for  elements of each column.

	int* temp= new int [n+1];//memory allocated for temporary vector array;
	
	
		for(int l=0;l<(n+1);l++)
		{	
			for(int m=0;m<(n+1);m++)
			{	
				temp[m]=0;

				for(int p=0;p<(n+1);p++)
					temp[m]=addTable(temp[m],mulTable(matrix_MT[p][l],matrix_MQ[p][m]));

			}

			for(int q=0;q<(n+1);q++)
			{
				matrix_MP[l][q]=0;

				for(int r=0;r<(n+1);r++)
					matrix_MP[l][q]=addTable(matrix_MP[l][q],mulTable(temp[r],matrix_MT[r][q]));
			}

		}

}

// convert MP into upper triangular matrix
void generate_matrix_MP_UT()
{
	for(int i=0;i<(n+1);i++)
		for(int j=i+1;j<(n+1);j++)
		{
			matrix_MP[i][j]=addTable(matrix_MP[i][j],matrix_MP[j][i]);
			matrix_MP[j][i]=0;
		}
}

// Macauley matrix of the public key
void generate_matrix_MPK()
{	
	matrix_MPK = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int j = 0 ; j < oilvariables ; j++ )
		matrix_MPK[j] = new int[((n+1)*(n+2))/2];//memory allocated for  elements of each column.
	
	generatepolyeqs();
	generate_matrix_MT();
	generate_L_inverse();
	
	matrix_MF = new int *[oilvariables] ;//memory allocated for elements of rows.
	

	for(int i=0;i<oilvariables;i++)
	{
		generate_matrix_MQ(i);
		generate_matrix_MF(i);
			
		generate_matrix_MP();
			
		generate_matrix_MP_UT();
		
		int counter=0;

		for(int k=0;k<n;k++)
			for(int l=k;l<n;l++)
			{
				matrix_MPK[i][counter]=matrix_MP[k][l];
				counter=counter+1;
			}

		for(int k=0;k<(n+1);k++)
		{
			matrix_MPK[i][counter]=matrix_MP[k][n];
			counter=counter+1;
		}
	}
	
	
}
// print and store output ======================================================================
void print_info()
{	
	myfile_pk.open("/home/ishtiaq/uov_output/Public Key.txt");
	myfile_pk<<"\t\t\t********Public Key Generation********\n\n";
	myfile_pk<<"Instruction\n\n";
	myfile_pk<<"1. If matrix has more than 15 elements in a column than colums are divided into blocks of 15 elements\n   and new elements start from next line.\n\n";
	myfile_pk<<"2. New row always start one line below the previous row.\n\n";
	myfile_pk<<"3. Matrix MT,Matrix MF,Matrix MPK and Public Key in Polynomial form are printed.\n\n\n";
	myfile_pk<<"Steps of Public Key Generation\n\n";
	myfile_pk<<"1. Matrix MT of size (n+1)*(n+1)is generated by using coefficients t and constants v of Linear Map T.\n\n";
	myfile_pk<<"2. Matrix MF is generated ny using secret coefficients (a,b,c and d) and constants(e)\n\n";
	myfile_pk<<"3. Matrix MQ of size (n+1)*(n+1) is generated by using coefficients(a,b,c and d) and constants(e) of Quadratic Polynomial.\n\n";
	myfile_pk<<"4. Matrix MP=Transpose(MT)*MQ*MT\n\n";
	myfile_pk<<"5. Matrix MQ and MP are generated oilvariables time.\n\n";
	myfile_pk<<"6. After generating matrix MP it is converted to triangular form\n\n";
	myfile_pk<<"7. Matrix MPK(Public Key) is generated using matrix MP(triangular form)\n\n";
	myfile_pk.close();

}
void print_matrix_MT()
{		
	myfile_pk.open("/home/ishtiaq/uov_output/Public Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MT is:\n\n";

	for(int i=0;i<(n+1);i++)
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
	myfile_pk.open("/home/ishtiaq/uov_output/Public Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MF  is:\n\n";

	for(int k=0;k<oilvariables;k++)
	{
		for(int j=0;j<((vinegarvariables*(n+oilvariables+1)+2*n+2)/2);j++)
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
	myfile_pk.open("/home/ishtiaq/uov_output/Public Key.txt",ios::app);
	myfile_pk<<"\n\nMatrix MPK is:\n\n";

	for(int k=0;k<oilvariables;k++)
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
void print_publickey()
{	
	print_info();
	print_matrix_MT();
	print_matrix_MF();
	print_matrix_MPK();


	myfile_pk.open("/home/ishtiaq/uov_output/Public Key.txt",ios::app);
	myfile_pk<<"\n\nPublic Key is written in polynomail form\n\n";
	
	for(int i=0;i<oilvariables;i++)
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
void print_privatekey()
{
	int l=0;
	
	myfile_pk.open("/home/ishtiaq/uov_output/Private Key.txt");
	
	myfile_pk<<"**********Private Key in Polynomail form***********\n";
	

	myfile_pk<<"\n\nQuadratic Polynomials are:\n";
	for(int i=0;i<oilvariables;i++)//Writing  Linerazied Polynomial equations=No of Oil Variables
	{   
		myfile_pk<<endl<<endl;
		myfile_pk<<"f("<<i+1<<"):\t";
			
		
		
		l=0;
		for(int j=0;j<vinegarvariables;j++)
                	for(int k=j;k<n;k++)
			{
				myfile_pk<<matrix_MF[i][l]<<"*x"<<j+1<<"*x"<<k+1<<"\t";
				if (matrix_MF[i][l]<10 && k<9 && j<9)
					myfile_pk<<"\t";
				
				myfile_pk<<"+\t";
				if((l+1)%5==0)
					myfile_pk<<endl<<"\t";
				l++;
			}

		
		for(int j=0;j<vinegarvariables;j++)//Writing coefficients c with vinegar variables
		{	
			myfile_pk<<coeffc[i][j]<<"*x"<<j+1<<"\t\t"<<"+"<<"\t";
			if((l+1)%5==0)
				myfile_pk<<endl<<"\t";			
			
			l++;
		}
		

		for(int j=0;j<oilvariables;j++)//Writing coefficients d with oil variables
		{	
			myfile_pk<<coeffd[i][j]<<"*x"<<j+1+vinegarvariables<<"\t\t"<<"+"<<"\t";
			if((l+1)%5==0)
				myfile_pk<<endl<<"\t";	
			
			l++;
			
		}
		myfile_pk<<coeffe[i];	//writing coefficients e
		
	
				
	}
	
	myfile_pk<<"\n\nLinear map T:\n";
		
	
	for(int i=0;i<oilvariables+vinegarvariables;i++)
	{	
		l=0;
		myfile_pk<<endl;
		myfile_pk<<"T("<<i+1<<"):\t";
		for (int j=0;j<oilvariables+vinegarvariables;j++)
		{
			myfile_pk<<coefft[i][j]<<"*w"<<j+1<<"\t+\t";
			if((l+1)%9==0)
				myfile_pk<<endl<<"\t";	
			l++;
		}
	
		myfile_pk<<coeffv[i]<<endl;
	}
	myfile_pk.close();
}
void store_privatekey()
{
	myfile_pk.open("/home/ishtiaq/data/Private Key UOV.txt");
	myfile_pk<<"o\n";
	myfile_pk<<oilvariables<<"\n";
	
	myfile_pk<<"v\n";
	myfile_pk<<vinegarvariables<<"\n";
	
	myfile_pk<<"MF\n";
	for(int k=0;k<oilvariables;k++)
	{
		for(int j=0;j<D;j++)
			myfile_pk<<matrix_MF[k][j]<<" ";
	        myfile_pk<<endl;
	}
	
	myfile_pk<<"c\n";
	for(int k=0;k<oilvariables;k++)
	{
		for(int j=0;j<vinegarvariables;j++)
			myfile_pk<<coeffc[k][j]<<" ";
	        myfile_pk<<endl;
	}
	
	myfile_pk<<"d\n";
	for(int k=0;k<oilvariables;k++)
	{
		for(int j=0;j<oilvariables;j++)
			myfile_pk<<coeffd[k][j]<<" ";
	        myfile_pk<<endl;
	}
	
	myfile_pk<<"e\n";
	for(int k=0;k<oilvariables;k++)
		myfile_pk<<coeffe[k]<<" ";
	        
	myfile_pk<<endl;
	
	myfile_pk<<"t\n";
	for(int k=0;k<n;k++)
	{
		for(int j=0;j<n;j++)
			myfile_pk<<coefft_inv[k][j]<<" ";
	        myfile_pk<<endl;
	}
	myfile_pk<<"v\n";
	for(int k=0;k<n;k++)
		myfile_pk<<coeffv[k]<<" ";
	        
	myfile_pk<<endl;
	
	myfile_pk.close();
			
		
	

	
}
void store_publickey()
{
	myfile_pk.open("/home/ishtiaq/data/Public Key UOV.txt");
	myfile_pk<<"o\n";
	myfile_pk<<oilvariables<<"\n";
	
	myfile_pk<<"v\n";
	myfile_pk<<vinegarvariables<<"\n";
	
	myfile_pk<<"MPK\n";
	for(int k=0;k<oilvariables;k++)
	{
		for(int j=0;j<((n+1)*(n+2))/2;j++)
			myfile_pk<<matrix_MPK[k][j]<<" ";
		myfile_pk<<endl;
	}
	myfile_pk.close();
	

}

