#include "loadkeysignhash.h"


int** matrix_MPK;

int* hashvalue;
int* message;
using namespace std;

void initialize(int oilvariables,int vinegarvariables)
{
	n=oilvariables+vinegarvariables;
	D = (vinegarvariables*(2*oilvariables+vinegarvariables+1))/2;
	
	
	matrix_MPK = new int *[oilvariables] ;//memory allocated for elements of rows.
	
	for( int j = 0 ; j < oilvariables ; j++ )
		matrix_MPK[j] = new int[((n+1)*(n+2))/2];//memory allocated for  elements of each column.
	
	hashvalue = new int [oilvariables] ;
	
	
	message = new int [oilvariables+vinegarvariables];//memory allocated for elements of vector.
	
}
void load_key_sign_hash()
{
  string line;
  ifstream myfile ("/home/ishtiaq/data/Public Key UOV.txt");
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
    	if (line.compare("MPK") == 0)
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
    				matrix_MPK[i][counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		}
	
	}
	else
		cout<<"Matrix MPK not present\n";

	    myfile.close();
  }

  else cout << "Unable to open Public Key file"; 

 
 ifstream myfile1 ("/home/ishtiaq/data/SignatureandHash_UOV.txt");
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
    				message[counter]=atoi(p);
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
