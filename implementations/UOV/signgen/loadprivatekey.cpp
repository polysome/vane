#include "loadprivatekey.h"


int** matrix_MF;
int** coeffc;
int** coeffd;
int* coeffe;
int** coefft;
int* coeffv;
int a;
using namespace std;

void initialize(int oilvariables,int vinegarvariables)
{
	n=oilvariables+vinegarvariables;
	D = (vinegarvariables*(2*oilvariables+vinegarvariables+1))/2;
	
	
	matrix_MF = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < (oilvariables) ; i++ )
		matrix_MF[i] = new int[D];//memory allocated for  elements of each column.
	
	coeffc = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < oilvariables ; i++ )
		coeffc[i] = new int[vinegarvariables];//memory allocated for  elements of each column.
 	
	coeffd = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < oilvariables ; i++ )
		coeffd[i] = new int[oilvariables];//memory allocated for  elements of each column.
	
	coeffe = new int [oilvariables] ;
	
	coefft = new int *[oilvariables+vinegarvariables] ;//memory allocated for elements of rows.
		
	for( int i = 0 ; i < oilvariables+vinegarvariables ; i++ )//memory allocated for  elements of each column.
		coefft[i] = new int[oilvariables+vinegarvariables];

	coeffv = new int [oilvariables+vinegarvariables];//memory allocated for elements of vector.
	
}
void load_privatekey()
{	
y: cout<<"Press 1 for generating signature by using private key for cyclic UOV,2 for LRS UOV and 3 for UOV: ";
  
  cin>>a; 
  string line;
  ifstream myfile;
  
 if(a==1)
	myfile.open("/home/ishtiaq/data/Private Key Cyclic.txt");
 if(a==2)
	myfile.open("/home/ishtiaq/data/Private Key LRS.txt");
 if(a==3)
	myfile.open("/home/ishtiaq/data/Private Key UOV.txt");
 if(a!=1 && a!=2 && a!=3)
 {
	cout<<"\nPress valid value\n";
	goto y;
 }


  if (myfile.is_open())
  {	
    	getline (myfile,line);
    	if (line.compare("o") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Oil Variables not present\n";
    	oilvariables= atoi(line.c_str());
	
	getline (myfile,line);
    	if (line.compare("v") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Vinegar Variables not present\n";
	vinegarvariables= atoi(line.c_str());
	
	initialize(oilvariables,vinegarvariables);
	
	getline (myfile,line);
    	if (line.compare("MF") == 0)
	{	for(int i=0;i<oilvariables;i++)
		{
			getline (myfile,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				matrix_MF[i][counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		}
	
	}
	else
		cout<<"Matrix MF not present\n";
	
	getline (myfile,line);
    	if (line.compare("c") == 0)
	{	for(int i=0;i<oilvariables;i++)
		{
			getline (myfile,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				coeffc[i][counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		}
	
	}
	else
		cout<<"Coeffc not present\n";
	
	getline (myfile,line);
    	if (line.compare("d") == 0)
	{	for(int i=0;i<oilvariables;i++)
		{
			getline (myfile,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				coeffd[i][counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		}
	
	}
	else
		cout<<"Coeffd not present\n";
	
	getline (myfile,line);
    	if (line.compare("e") == 0)
	{	
			getline (myfile,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				coeffe[counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		
	
	}
	else
		cout<<"Coeffe not present\n";

	getline (myfile,line);
    	if (line.compare("t") == 0)
	{	for(int i=0;i<n;i++)
		{
			getline (myfile,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				coefft[i][counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		}
	
	}
	else
		cout<<"Coefft not present\n";

	getline (myfile,line);
    	if (line.compare("v") == 0)
	{	
			getline (myfile,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				coeffv[counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		
	
	}
	else
		cout<<"Coeffv not present\n";
	

	
    
    myfile.close();
  }

  else cout << "Unable to open file"; 
}
