#include "signvericyclic.h"

using namespace std;

ofstream myfile_sv;


int* hash_verify=0;
int result;



void signverificationcyclic()
{	
	int counter=0;
	int minor=v2-1;

	hash_verify=new int[m];

	int* temp=new int[n+1];
	
	int** a = new int *[2];	
		for( int i = 0 ; i < 2 ; i++ )
			a[i] = new int[n+1];
	
	for (int j=0;j<v1;j++)
	{
		a[0][j]=0;
		
		for(int k=0;k<=j;k++)
			a[0][j]=addTable(a[0][j],mulTable(signature[k],matrix_PK[0][k][j]));

		temp[j]=a[0][j];
		
	}

	for(int j=v1;j<minor;j++)
	{
		a[0][j]=0;
		
		for(int k=0;k<v1;k++)
			a[0][j]=addTable(a[0][j],mulTable(signature[k],matrix_PK[0][k][j]));

		temp[j]=a[0][j];

		for(int k=v1;k<=j;k++)
			temp[j]=addTable(temp[j],mulTable(signature[k],matrix_PK[0][k][j]));
	}

	for(int j=minor;j<n;j++)
	{
		temp[j]=0;
		
		for(int k=0;k<=j;k++)
			temp[j]=addTable(temp[j],mulTable(signature[k],matrix_PK[0][k][j]));		
	
	}	
	
	temp[n]=0;
	for(int k=0;k<n;k++)
		temp[n]=addTable(temp[n],mulTable(signature[k],matrix_PK[0][k][n]));

	temp[n]=addTable(temp[n],matrix_PK[0][n][n]);
	
	hash_verify[0]=0;
	
	for(int j=0;j<n;j++)
			hash_verify[0]=addTable(hash_verify[0],mulTable(signature[j],temp[j]));
		
	hash_verify[0]=addTable(hash_verify[0],temp[n]);

	if(hash_verify[0]==hashvalue[0])
		counter++;
	else
		goto q;


	for(int l=1;l<o1;l++)
	{	
		
		for(int j=v2;j<n;j++)
		{
			temp[j]=0;
		
			for(int k=0;k<=j;k++)
				temp[j]=addTable(temp[j],mulTable(signature[k],matrix_PK[l][k][j]));	
		}
		temp[n]=0;
		for(int k=0;k<n;k++)
			temp[n]=addTable(temp[n],mulTable(signature[k],matrix_PK[l][k][n]));

		temp[n]=addTable(temp[n],matrix_PK[l][n][n]);

		temp[minor]=a[0][minor-1];
		for(int k=v1;k<v2;k++)
			temp[minor]=addTable(temp[minor],mulTable(signature[k],matrix_PK[l][k][minor]));

		for(int j=v2-2;j>=v1;j--)
		{
			a[0][j]=a[0][j-1];
			temp[j]=a[0][j];
			
			for(int k=v1;k<=j;k++)
				temp[j]=addTable(temp[j],mulTable(signature[k],matrix_PK[l][k][j]));
			
		
		}
	
		for(int j=v1-1;j>=1;j--)
		{
			a[0][j]=addTable(a[0][j-1],mulTable(signature[j],matrix_PK[l][j][j]));
			temp[j]=a[0][j];
		}
	
		a[0][0]=mulTable(signature[0],matrix_PK[l][0][0]);
		temp[0]=a[0][0];

		hash_verify[l]=0;
	
		for(int j=0;j<n;j++)
			hash_verify[l]=addTable(hash_verify[l],mulTable(signature[j],temp[j]));
		
		hash_verify[l]=addTable(hash_verify[l],temp[n]);

		if(hash_verify[l]==hashvalue[l])
			counter++;
		else
			goto q;

		
	}
	
	temp[n-1]=0;
	for(int k=0;k<n;k++)
		temp[n-1]=addTable(temp[n-1],mulTable(signature[k],matrix_PK[o1][k][n-1]));
		
	
	temp[n]=0;
	for(int k=0;k<n;k++)
		temp[n]=addTable(temp[n],mulTable(signature[k],matrix_PK[o1][k][n]));

	temp[n]=addTable(temp[n],matrix_PK[o1][n][n]);


	for(int j=v2;j<n-1;j++)
	{
		a[1][j-v2+1]=0;
		for(int k=0;k<v2;k++)
			a[1][j-v2+1]=addTable(a[1][j-v2+1],mulTable(signature[k],matrix_PK[o1][k][j]));

		temp[j]=a[1][j-v2+1];

		for(int k=v2;k<=j;k++)
			temp[j]=addTable(temp[j],mulTable(signature[k],matrix_PK[o1][k][j]));
		

	}

	a[1][0]=0;
	for(int k=v1;k<v2;k++)
		a[1][0]=addTable(a[1][0],mulTable(signature[k],matrix_PK[o1][k][minor]));

	temp[minor]=addTable(a[0][minor-1],a[1][0]);	

	for(int j=v2-2;j>=v1;j--)
	{
		a[0][j]=a[0][j-1];
		for(int k=v1;k<=j;k++)
			a[0][j]=addTable(a[0][j],mulTable(signature[k],matrix_PK[o1][k][j]));	
		
		temp[j]=a[0][j];
	}

	for(int j=v1-1;j>=1;j--)
	{
		a[0][j]=addTable(a[0][j-1],mulTable(signature[j],matrix_PK[o1][j][j]));
		temp[j]=a[0][j];	
	}

	a[0][0]=mulTable(signature[0],matrix_PK[o1][0][0]);
	temp[0]=a[0][0];
	
	hash_verify[o1]=0;
	
	for(int j=0;j<n;j++)
		hash_verify[o1]=addTable(hash_verify[o1],mulTable(signature[j],temp[j]));
	
	hash_verify[o1]=addTable(hash_verify[o1],temp[n]);

	if(hash_verify[o1]==hashvalue[o1])
		counter++;
	else
		goto q;


	for(int l=o1+1;l<m;l++)
	{

		temp[n]=0;
		for(int k=0;k<n;k++)
			temp[n]=addTable(temp[n],mulTable(signature[k],matrix_PK[l][k][n]));

		temp[n]=addTable(temp[n],matrix_PK[l][n][n]);

		temp[n-1]=a[1][o2-1];
		for(int k=v2;k<n;k++)
			temp[n-1]=addTable(temp[n-1],mulTable(signature[k],matrix_PK[l][k][n-1]));

		
		for(int j=n-2;j>=v2+1;j--)
		{
			a[1][j-v2+1]=a[1][j-v2];
			temp[j]=a[1][j-v2+1];

			for(int k=v2;k<=j;k++)
				temp[j]=addTable(temp[j],mulTable(signature[k],matrix_PK[l][k][j]));
			

		}

		a[1][1]=a[1][0];
		for(int k=0;k<v1;k++)
			a[1][1]=addTable(a[1][1],mulTable(signature[k],matrix_PK[l][k][v2]));

		temp[v2]=addTable(a[1][1],mulTable(signature[v2],matrix_PK[l][v2][v2]));

		a[1][0]=0;
		for(int k=v1;k<v2;k++)
			a[1][0]=addTable(a[1][0],mulTable(signature[k],matrix_PK[l][k][minor]));

		temp[minor]=addTable(a[0][minor-1],mulTable(signature[minor],matrix_PK[l][minor][minor]));
	
		for(int j=v2-2;j>=1;j--)
		{
			a[0][j]=addTable(a[0][j-1],mulTable(signature[j],matrix_PK[l][j][j]));
			temp[j]=a[0][j];
		}
		
		a[0][0]=mulTable(signature[0],matrix_PK[l][0][0]);
		temp[0]=a[0][0];

		hash_verify[l]=0;
	
		for(int j=0;j<n;j++)
			hash_verify[l]=addTable(hash_verify[l],mulTable(signature[j],temp[j]));
		
		hash_verify[l]=addTable(hash_verify[l],temp[n]);

		if(hash_verify[l]==hashvalue[l])
			counter++;
		else
			goto q;

	}


	
q:	if(counter==m)
		result=1;
	else
		result=0;
	
	
}

void print_signatureverification(int result)
{	
	myfile_sv.open("/home/ishtiaq/rainbow/rainbow_cyclic_output/Signature Verification.txt");
	myfile_sv<<"\n\nSignature Verification\n\n";
	myfile_sv<<"Hash values generated by verification process\n\n";
	myfile_sv<<"Hashvalue= (";
	for(int i=0;i<m;i++)
	{	
		myfile_sv<<hash_verify[i];
		if(i!=(m-1))
			myfile_sv<<", ";
	}
	myfile_sv<<")"<<endl;

	if(result==1)
	{
		myfile_sv<<"\n\nHash values generated by verification process are equal to original Hashvalues\nSo Signature Authenticated";
	}
	else
	{
		myfile_sv<<"\n\nHash values generated by verification process are not equal to original Hashvalues\nSo Signature Authentication Failed";
	}
	myfile_sv.close();

}
