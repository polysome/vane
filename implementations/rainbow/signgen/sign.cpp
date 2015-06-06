#include "basic.h"
#include "sign.h"



 int* hashvalue;
 int* solution_linearmapS;
 int* vinegarVariable;
 int* solution;//It stores oil and vinegar variables
 int* signature;//It contains signature of message

 int** linearmap_S;
 int** linearmap_T;
 int* temp;


 using namespace std;
 ofstream myfile_sg;

void randomvalues4vinegarVariables()
{	
	vinegarVariable = new int [v1] ;//memory allocated for elements of rows.
	
	for(int index=0; index<v1; index++)
		vinegarVariable[index] = (rand()%256); 
		
}

void hashvaluegenerator()
{
	hashvalue = new int [m] ;//memory allocated for elements of rows.
	randvectorgenerator(m,hashvalue);
}
void linearmap_of_S()
{	
	temp=new int[m];
	solution_linearmapS=new int [m];

	for(int i=0;i<m;i++)
		temp[i]=addTable(matrix_MS_inverse[i][m],hashvalue[i]);

	for(int i=0;i<m;i++)
	{	
		solution_linearmapS[i]=0;
		for(int j=0;j<m;j++)
			solution_linearmapS[i]=addTable(solution_linearmapS[i],mulTable(matrix_MS_inverse[i][j],temp[j]));
	}


}


void Quadtraticmap()
{	
	int counter,c_start=(v1*(v1+1))/2+o1*v1,d_start=c_start+v1,e_start=d_start+o1;
	int counter2,c_start2=(v2*(v2+1))/2+o2*v2,d_start2=c_start2+v2,e_start2=d_start2+o2;
	
	int **polyeq_1;
	polyeq_1 = new int *[o1] ;//memory allocated for elements of rows.
	for( int i = 0 ; i < o1 ; i++ )
		polyeq_1[i] = new int[o1+1];//memory allocated for  elements of each column.

	int **polyeq_2;
	polyeq_2 = new int *[o2] ;//memory allocated for elements of rows.
	for( int i = 0 ; i < o2 ; i++ )
		polyeq_2[i] = new int[o2+1];//memory allocated for  elements of each column.

	solution = new int [n] ;//memory allocated for elements of rows.
	
	
a:	randomvalues4vinegarVariables();// random value generation function
	

	for(int i=0;i<o1;i++)
	{	
		for(int j=0;j<o1;j++)
		{	
			counter=v1;
				
			polyeq_1[i][j]=matrix_MF[i][d_start+j];
			
			for(int k=0;k<v1;k++)
			{
				polyeq_1[i][j]=addTable(polyeq_1[i][j],mulTable(matrix_MF[i][j + counter ],vinegarVariable[k]));
				counter=counter+v2-k-1;
			}   
		
		}		


		counter=0;
		polyeq_1[i][o1]=matrix_MF[i][e_start];

		
		for(int j=0;j<v1;j++)
		{	
			polyeq_1[i][o1]=addTable(polyeq_1[i][o1],mulTable(matrix_MF[i][c_start+j],vinegarVariable[j]));

			
			for(int k=j;k<v1;k++)
			{	
				int temp2=mulTable(vinegarVariable[j],vinegarVariable[k]);
				polyeq_1[i][o1]=addTable(polyeq_1[i][o1],mulTable(matrix_MF[i][counter],temp2));
				counter=counter+1;
			}
			counter=counter+o1;
		}
	}
	
	for(int i=0;i<o1;i++)//Adding solution of Linear Map S to constants
		polyeq_1[i][o1]=addTable(polyeq_1[i][o1],solution_linearmapS[i]);
	
		
	//Calling function of Guass elimination to solve polynomial equations
	int* guasspointer=guasselimination(polyeq_1,o1,o1+1);
	
	if(guasspointer[0]==-1)//If solution doesnt exist choose new values for vinegar variables			
		goto a;
		
	for(int i=0;i<v1;i++)//Assigning value of vinegar variables to solution array
		solution[i]=vinegarVariable[i];
			
	for(int j=0;j<o1;j++)//Assigning value of oil variables to solution array
		solution[v1+j]=guasspointer[j];
	
		
	for(int i=0;i<o2;i++)
	{	
		for(int j=0;j<o2;j++)
		{	
			counter2=v2;
				
			polyeq_2[i][j]=matrix_MF[o1+i][d_start2+j];
			
			for(int k=0;k<v2;k++)
			{
				polyeq_2[i][j]=addTable(polyeq_2[i][j],mulTable(matrix_MF[o1+i][j + counter2 ],solution[k]));
				counter2=counter2+n-k-1;
			}   
		
		}		


		counter2=0;
		polyeq_2[i][o2]=matrix_MF[o1+i][e_start2];

		
		for(int j=0;j<v2;j++)
		{	
			polyeq_2[i][o2]=addTable(polyeq_2[i][o2],mulTable(matrix_MF[o1+i][c_start2+j],solution[j]));

			
			for(int k=j;k<v2;k++)
			{	
				int temp2=mulTable(solution[j],solution[k]);
				polyeq_2[i][o2]=addTable(polyeq_2[i][o2],mulTable(matrix_MF[o1+i][counter2],temp2));
				counter2=counter2+1;
			}
			counter2=counter2+o2;
		}
	}
	for(int i=0;i<o2;i++)//Adding solution of Linear Map S to constants
		polyeq_2[i][o2]=addTable(polyeq_2[i][o2],solution_linearmapS[o1+i]);
	
	

	int* guass_pointer=guasselimination(polyeq_2,o2,o2+1);
	
	if(guass_pointer[0]==-1)//If solution doesnt exist choose new values for vinegar variables			
		goto a;
			
	for(int j=0;j<o2;j++)//Assigning value of oil variables to solution array
		solution[v2+j]=guass_pointer[j];
	
		
}
void linearmap_of_T()
{	
	
	temp=new int[n];
	signature=new int [n];

	for(int i=0;i<n;i++)
		temp[i]=addTable(matrix_MT_inverse[i][n],solution[i]);

	for(int i=0;i<n;i++)
	{	
		signature[i]=0;
		for(int j=0;j<n;j++)
			signature[i]=addTable(signature[i],mulTable(matrix_MT_inverse[i][j],temp[j]));
	}
	

}

int* guasselimination(int** ge,int rows,int columns)
{	
	int *solnarray = 0;
	solnarray = new int [rows] ;//memory allocated for elements of rows.

	int phajja;
	int maxy=rows-1;
	int value=0;
	int raxy=rows-1;
	
	int* divisionPointer = divisionlookuparray();
	int counter=0;
	
	for(int l=0;l<rows-1;l++)//reducing the matrix to upper triangular form using elemntary row operations,also checking necssary conditions at which solution doesn't exist.
	{	
		while(ge[counter][l]==0&&counter<rows)
		{   
			counter++;
		}
			
		if(counter==rows)
		{
			phajja=0;
			goto m;
		}
		if(counter!=rows)
		{
			for(int c=0;c<columns;c++)
			{
				int tempelement=ge[l][c];
				ge[l][c]=ge[counter][c];
				ge[counter][c]=tempelement;

			}
			counter=0;
		}

		counter=l+1;
		
		for(int k=0;k<raxy;k++)
		{
			for(int i=0;i<256;i++)
			{
				int result=addTable(ge[l+k+1][l],mulTable(ge[l][l],i));
				if(result==0)
				{
					value=i;
					break;
				}	
			}
	
			for(int i=0;i<columns;i++)
			{
				int temp=mulTable(ge[l][i],value);
				ge[l+k+1][i]=addTable(ge[l+k+1][i],temp);
			}
	
			int abc=0;
			while(ge[l+k+1][abc]==0&&abc<rows)
			{
				abc++;
			}
			if(abc==rows)
			{
				phajja=0;
				goto m;
							
			}
				
		}
		raxy--;
	}
	
	

	for(int l=rows-1;l>0;l--)//reducing the upper triangular matrix to diagonal form by using elementary row oprations,also checking conditions at which soln doesn't exist.
	{
		for(int k=0;k<maxy;k++)
		{
			for(int i=0;i<256;i++)
			{
				int result=addTable(ge[l-k-1][l],mulTable(ge[l][l],i));
				if(result==0)
				{
					value=i;
					break;
				}	
			}
			for(int i=0;i<columns;i++)
			{
				int temp=mulTable(ge[l][i],value);
				ge[l-k-1][i]=addTable(ge[l-k-1][i],temp);
			}
		}
				
		maxy--;
	}
	
	
	for(int i=0;i<rows;i++)//dividing to make the diagonal of matrix to one
	{
		int temp=ge[i][i];
		for(int j=0;j<columns;j++)
		{
			ge[i][j]=divideTable(ge[i][j],temp,divisionPointer);
		}
		solnarray[i]=ge[i][columns-1];
		
	}
	
	phajja=1;
	
m: if(phajja==0)
	{
		solnarray[0]=-1;
	}
	
	return solnarray;
}



void generatesignature()
{
	hashvaluegenerator();
	linearmap_of_S();
	Quadtraticmap();
	linearmap_of_T();
}

void print_sign_generation()
{	
	//int l=0;
	if(a==1)
	myfile_sg.open("/home/ishtiaq/rainbow/rainbow_conventional_output/Signature Generation.txt");
	if(a==2)
	myfile_sg.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Signature Generation.txt");
	if(a==3)
	myfile_sg.open("/home/ishtiaq/rainbow/rainbow_LRS_output/Signature Generation.txt");
	
	
	
	myfile_sg<<"**********Signature Generation***********\n";
	
	
	myfile_sg<<"\n\nHashvalue= (";
	for(int i=0;i<m;i++)
	{	
		myfile_sg<<hashvalue[i];
		if(i!=(m-1))
			myfile_sg<<", ";
	}
	myfile_sg<<")"<<endl;
	
	myfile_sg<<"\n\nVinegar Variables"<<endl<<endl;
	for(int i=0;i<v1;i++)
		myfile_sg<<"x"<<i+1<<"="<<solution[i]<<"\t\t";
	
	myfile_sg<<"\n\nOil Variables(O1)"<<endl<<endl;
	for(int i=v1;i<v2;i++)
		myfile_sg<<"x"<<i+1<<"="<<solution[i]<<"\t\t";

	myfile_sg<<"\n\nOil Variables(O2)"<<endl<<endl;
	for(int i=v2;i<n;i++)
		myfile_sg<<"x"<<i+1<<"="<<solution[i]<<"\t\t";	

	myfile_sg<<"\n\nSignature = (";
	for(int i=0;i<n;i++)
		{	
			myfile_sg<<signature[i]; 
			if(i!=((n)-1))
				myfile_sg<<", ";
		}
	myfile_sg<<")"<<endl;
	myfile_sg.close();
	
}

void store_hash_sign()
{
	if(a==1)
	myfile_sg.open("/home/ishtiaq/rainbow/data/SignatureandHash_conventional.txt");
	if(a==2)
	myfile_sg.open("/home/ishtiaq/rainbow/data/SignatureandHash_cyclic.txt");
	if(a==3)
	myfile_sg.open("/home/ishtiaq/rainbow/data/SignatureandHash_LRS.txt");
	
	myfile_sg<<"Hashvalue\n";
	for(int i=0;i<m;i++)
		myfile_sg<<hashvalue[i]<<" ";
	myfile_sg<<endl;
	myfile_sg<<"Signature\n";
	for(int i=0;i<n;i++)
		myfile_sg<<signature[i]<<" ";
	myfile_sg<<endl;
		
}

