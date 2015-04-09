#include "basic.h"
#include "sign.h"



 int* hashvalue;
 int* temp;
 int *solution;//It stores oil and vinegar variables
 int *message;//It contains signature of message


 using namespace std;
 ofstream myfile_sg;

int* randomvalues4vinegarVariables()
{	
	int *vinegarVariable = 0;
	vinegarVariable = new int [vinegarvariables] ;//memory allocated for elements of rows.
	
	for(int index=0; index<vinegarvariables; index++)
		vinegarVariable[index] = (rand()%255); 
		
	return vinegarVariable;
}






void hashvaluegenerator()
{
	hashvalue = new int [oilvariables] ;//memory allocated for elements of rows.
	randvectorgenerator(oilvariables,hashvalue);
}

int* solvepolynomialequations()
{	
	
	hashvaluegenerator();
	
	int **polyeq = 0;//It stores polynomial eq after linearized
	polyeq = new int *[oilvariables] ;//memory allocated for elements of rows.
	for( int i = 0 ; i < oilvariables ; i++ )
		polyeq[i] = new int[oilvariables+1];//memory allocated for  elements of each column.

	
	solution = new int [oilvariables+vinegarvariables] ;//memory allocated for elements of rows.
	
	
a:	int* randomPointer = randomvalues4vinegarVariables();// random value generation function
	

	for(int k=0;k<oilvariables;k++)//After putting value of vinegar variables adding coefficients of oil variables
		for(int j=0;j<oilvariables;j++)
		{
			polyeq[k][j]=coeffd[k][j];
			for(int i=0;i<vinegarvariables;i++)
			   polyeq[k][j]=addTable(polyeq[k][j],mulTable(matrix_MF[k][vinegarvariables + n*i - ((i*(i+1))/2)+ j ],randomPointer[i]));
		}
//
	    
	for(int i=0;i<oilvariables;i++)//After putting value of vinegar variables, adding constants
	{
		polyeq[i][oilvariables]=coeffe[i];
		
		for(int j=0;j<vinegarvariables;j++)
		{	
			polyeq[i][oilvariables]=addTable(polyeq[i][oilvariables],mulTable(coeffc[i][j],randomPointer[j]));
			
			for(int k=0;k<(vinegarvariables-j);k++)
			{	
				int temp2=mulTable(randomPointer[j],randomPointer[k+j]);
				polyeq[i][oilvariables]=addTable(polyeq[i][oilvariables],mulTable(matrix_MF[i][n*j-((j*(j-1))/2)+k],temp2));
			}
		}
	 }

	
	for(int i=0;i<oilvariables;i++)//Adding Hashvalues to constants
		polyeq[i][oilvariables]=addTable(polyeq[i][oilvariables],hashvalue[i]);
			
	int* guasspointer=guasselimination(polyeq,oilvariables,oilvariables+1);//Calling function of Guass elimination to solve polynomial equations
	
	if(guasspointer[0]==-1)//If solution doesnt exist choose new values for vinegar variables			
		goto a;
		
	for(int i=0;i<vinegarvariables;i++)//Assigning value of vinegar variables to solution array
		solution[i]=randomPointer[i];
			
	for(int j=0;j<oilvariables;j++)//Assigning value of oil variables to solution array
		solution[j+vinegarvariables]=guasspointer[j];
			
	return solution;
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


int* linearmapT(int* solution)
{	
		
	temp = new int [n] ;//memory allocated for elements of rows.
	

	
	message = new int [oilvariables+vinegarvariables];//memory allocated for elements of rows.
	for(int i=0;i<oilvariables+vinegarvariables;i++)//Adding v and oil and vinegar variables.
		temp[i]=addTable(coeffv[i],solution[i]);
		
	for(int i=0;i<n;i++)
	{	
		message[i]=0;
		for(int j=0;j<n;j++)
			message[i]=addTable(message[i],mulTable(coefft[i][j],temp[j]));
	
	}

	
	return message;

}



void generatesignature()
{
	int* solutionpointer=solvepolynomialequations();//After solving polynomial equations it returns vinegar and oil variables
	int* signmessage=linearmapT(solutionpointer);//After solving Linear Map T it signature of message
}

void print_sign_generation()
{	
	//int l=0;
    myfile_sg.open("/home/ishtiaq/uov_output/Signature Generation.txt");
	myfile_sg<<"**********Signature Generation***********\n";
	
	
	myfile_sg<<"\n\nHashvalue= (";
	for(int i=0;i<oilvariables;i++)
	{	
		myfile_sg<<hashvalue[i];
		if(i!=(oilvariables-1))
			myfile_sg<<", ";
	}
	myfile_sg<<")"<<endl;
	
	myfile_sg<<"\n\n1.Set Quadratic Polynomials f(i)=h(i)  (i=1,...oilvariables)\n";
	myfile_sg<<"2.Random values are assigned to the vinegar variables the polynomials f(i) become linear.\n";
	myfile_sg<<"3.Linear system is solved by Guass Elimination to get oilvariables,if solution doesn't exist\n  new random values assigned to vinegar variables.\n";
	
	myfile_sg<<"\n\nVinegar Variables"<<endl<<endl;
	for(int i=0;i<vinegarvariables;i++)
		myfile_sg<<"x"<<i+1<<"="<<solution[i]<<"\t\t";
	
	myfile_sg<<"\n\nOil Variables"<<endl<<endl;
	for(int i=vinegarvariables;i<oilvariables+vinegarvariables;i++)
		myfile_sg<<"x"<<i+1<<"="<<solution[i]<<"\t\t";
	
	myfile_sg<<"\n\n\n4.Set Linear Map T(i)=x(i) (i=1,...,Oilvariables+Vinegarvariables)";
	

	myfile_sg<<"\n\nSignature = (";
	for(int i=0;i<oilvariables+vinegarvariables;i++)
		{	
			myfile_sg<<message[i]; 
			if(i!=((oilvariables+vinegarvariables)-1))
				myfile_sg<<", ";
		}
	myfile_sg<<")"<<endl;
	myfile_sg.close();
	
}
void store_hash_sign()
{
    myfile_sg.open("/home/ishtiaq/data/SignatureandHash_UOV.txt");
	myfile_sg<<"Hashvalue\n";
	for(int i=0;i<oilvariables;i++)
		myfile_sg<<hashvalue[i]<<" ";
	myfile_sg<<endl;
	myfile_sg<<"Signature\n";
	for(int i=0;i<n;i++)
		myfile_sg<<message[i]<<" ";
	myfile_sg<<endl;
		
}
