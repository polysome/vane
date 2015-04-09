#include "loadprivatekey.h"

int a;
int** matrix_MF;
int** matrix_MS_inverse;
int** matrix_MT_inverse;

using namespace std;

void initialize(int o1,int o2,int v1,int v2,int m,int n)
{	
	matrix_MF = new int *[m] ;//memory allocated for elements of rows.
	
	for( int i = 0 ; i < o1 ; i++ )
	matrix_MF[i] = new int[v1*o1 + (v1*(v1+1))/2 + v2 + 1];//memory allocated for  elements of each column.
	
	for( int i = o1 ; i < m ; i++ )
	matrix_MF[i] = new int[v2*o2 + (v2*(v2+1))/2 + v2 + o2 + 1];//memory allocated for  elements of each column.
	
		
	matrix_MS_inverse = new int *[m] ;//memory allocated for elements of rows.
		
	for( int i = 0 ; i < m ; i++ )//memory allocated for  elements of each column.
		matrix_MS_inverse[i] = new int[(m+1)];

	matrix_MT_inverse = new int *[n] ;//memory allocated for elements of rows.
		
	for( int i = 0 ; i < n ; i++ )//memory allocated for  elements of each column.
		matrix_MT_inverse[i] = new int[(n+1)];
	
}
void load_privatekey()
{

 y: cout<<"Press 1 for generating signature by using private key for conventional Rainbow ,2 for cyclic Rainbow and 3 for LRS Rainbow: ";
  
  cin>>a; 
  string line;
  ifstream myfile;
  
 if(a==1)
	myfile.open("/home/ishtiaq/rainbow/data/PrivateKey.txt");
 if(a==2)
	myfile.open("/home/ishtiaq/rainbow/data/PrivateKeycyclic.txt");
 if(a==3)
	myfile.open("/home/ishtiaq/rainbow/data/PrivateKeyLRS.txt");
 if(a!=1 && a!=2 && a!=3)
 {
	cout<<"\nPress valid value\n";
	goto y;
 }	
 
 if (myfile.is_open())
  {	
    	getline (myfile,line);
    	if (line.compare("o1") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Oil Variables not present\n";
    	o1= atoi(line.c_str());

	getline (myfile,line);
    	if (line.compare("o2") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Oil Variables not present\n";
    	o2= atoi(line.c_str());

	getline (myfile,line);
    	if (line.compare("v1") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Oil Variables not present\n";
    	v1= atoi(line.c_str());

	getline (myfile,line);
    	if (line.compare("v2") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Oil Variables not present\n";
    	v2= atoi(line.c_str());

	getline (myfile,line);
    	if (line.compare("m") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Oil Variables not present\n";
    	m= atoi(line.c_str());

	getline (myfile,line);
    	if (line.compare("n") == 0)
		getline (myfile,line);
	else
		cout<<"Error Number of Oil Variables not present\n";
    	n= atoi(line.c_str());
	
	
	initialize(o1,o2,v1,v2,m,n);
	
	getline (myfile,line);
    	if (line.compare("MF") == 0)
	{	for(int i=0;i<m;i++)
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
    	if (line.compare("S") == 0)
	{	for(int i=0;i<m;i++)
		{
			getline (myfile,line);
			char * cstr, *p;
			int counter=0;
			cstr = new char [line.size()+1];
  			strcpy (cstr, line.c_str());
  			p=strtok (cstr," ");
  			while (p!=NULL)
  			{
    				matrix_MS_inverse[i][counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		}
	
	}
	else
		cout<<"Matrix MS inverse not present\n";

	getline (myfile,line);
    	if (line.compare("T") == 0)
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
    				matrix_MT_inverse[i][counter]=atoi(p);
				counter=counter+1;
    				p=strtok(NULL," ");
  			}

  			delete[] cstr; 
			
		}
	
	}
	else
		cout<<"Matrix MT inverse not present\n";


    myfile.close();
  }

  else cout << "Unable to open file"; 
}
