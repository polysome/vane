#include "loadkeysignhash.h"


int*** matrix_PK;

int* hashvalue;
int* signature;
using namespace std;

void initialize(int m,int n)
{

	matrix_PK = new int **[m];
  	for (int i = 0; i < m; i++) {
    		matrix_PK[i] = new int*[n+1];
		for (int j = 0; j < n+1; j++)
      			matrix_PK[i][j] = new int[n+1];
	}
	
	hashvalue = new int [m] ;
	
	
	signature = new int [n];//memory allocated for elements of vector.
	
}
void load_key_sign_hash()
{
  string line;
  ifstream myfile ("/home/ishtiaq/rainbow/data/PublicKey.txt");
  if (myfile.is_open())
  {
    	getline (myfile,line);
    	if (line.compare("m") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Oil Variables(o1+o2) not present\n";
    	m= atoi(line.c_str());
	
	getline (myfile,line);
    	if (line.compare("n") == 0)
		getline (myfile,line);
	else
		cout<<"Error value of n not present\n";
	n= atoi(line.c_str());
	
	initialize(m,n);
	
	
	getline (myfile,line);
    	if (line.compare("pk") == 0)
	{	
		for(int k=0;k<m;k++)
		{		
		for(int i=0;i<n+1;i++)
		{
			getline (myfile,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				matrix_PK[k][i][counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		}
		}
	
	}
	else
		cout<<"Matrix PK not present\n";

	    myfile.close();
  }

  else cout << "Unable to open Public Key file"; 

 
 ifstream myfile1 ("/home/ishtiaq/rainbow/data/SignatureandHash_conventional.txt");
  if (myfile1.is_open())
  {
	
	getline (myfile1,line);
    	if (line.compare("Hashvalue") == 0)
	{	
			getline (myfile1,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				hashvalue[counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		
	
	}
	else
		cout<<"Hash Values not present\n";

	getline (myfile1,line);
    	if (line.compare("Signature") == 0)
	{	
			getline (myfile1,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				signature[counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		
	
	}
	else
		cout<<"Signature not present\n";

	
    
    myfile1.close();
  }

  else cout << "Unable to open file"; 
}
