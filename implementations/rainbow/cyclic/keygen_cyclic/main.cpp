// Ishtiaqcourage.cpp : Defines the entry point for the console application.
//

#include "basic.h"
#include "keygen.h"
#include "Timer1.h"
#include<sys/stat.h>
#include<sys/types.h>



using namespace std;



int main() 
{	

	Timer1 timer;
	char c;
	const int length=100;
	double time_measure[length];
	double *result_time=0;
        unsigned long long cpu_cycles[length];
        unsigned long long* result_cpu;
        unsigned long long cycles;
	

	cout<<"\n**********Implementation of Key Generation of Cyclic Rainbow Scheme***********\n";
	
x:	cout<<"\nPlease enter number of V1 (Vinegar Variables e-g 17)  ";
	cin>>v1;	
	cout<<"\nPlease enter number of O1 (Oil Variables e-g 13)  ";
	cin>>o1;
	cout<<"\nPlease enter number of O2 (Oil Variables e-g 13)  ";
	cin>>o2;
	
     	cout<<"\nIf you want change input press 'c',If you want to exit press 'e',\nelse press any alphabet to continue  ";
    	cin>>c;
	if(c=='c')
		goto x;
	if(c=='e')
		goto y;
	
	m=o1+o2;
	n=v1+o1+o2;
        v2=v1+o1;
	D = ((n+1)*(n+2))/2,D1=(v1*(v1+1))/2 + v1*o1,D2=(v2*(v2+1))/2 + v2*o2;
        D21=D2-D1;
	
	mkdir("/home/ishtiaq/rainbow",0777);
	mkdir("/home/ishtiaq/rainbow/rainbow_cyclic_output",0777);
	mkdir("/home/ishtiaq/rainbow/data",0777);
	srand((unsigned)time(0));
	
    	cout<<"\nGenerating Key ...\n";
	timer.start();// Start the timer
	cycles = rdtsc();
	generate_key();
	cycles = rdtsc() - cycles;  
	timer.stop();// Stop the timer
	cout<<"Generated\n";
	cout<<"Time Elapsed in Seconds for Generating Public Key: "<<timer.getElapsedTimeInSec(); 
	cout<<"\nNumber of Cpu Cycles used for Generating Public Key: "<<cycles;
	
	cout<<"\nWriting the matrix B1 and B2 generated during key generation\nPrivate Key in Matrix(MS,MT and MF) and Public Key in both \nMatrix(MPK) and Polynomial form in home\\ishtiaq\\rainbow\\rainbow_cyclic_output\\Key.txt...\n";
	print_key();//It writes data in text file
	cout<<"Completed\n";

	store_publickey();
	store_privatekey();
	
	
y:	return 0;
}




 


