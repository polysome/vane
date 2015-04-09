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
        
	
	cout<<"\nImplementation of Unbalanced Oil and Vinegar Scheme\n";
	
x:	cout<<"\nPlease enter number of Oil Variables(e-g 2)  ";
	cin>>oilvariables;
	cout<<"\nPlease enter number of Vinegar Variables(e-g 4)  ";
	cin>>vinegarvariables;
        cout<<"\nIf you want change input press 'c',If you want to exit press 'e',\nelse press any alphabet to continue  ";
        cin>>c;
		if(c=='c')
			goto x;
		if(c=='e')
			goto y;


	n=oilvariables+vinegarvariables;
	D = (vinegarvariables*(2*oilvariables+vinegarvariables+1))/2;

	mkdir("/home/ishtiaq/uov_output",0777);
	mkdir("/home/ishtiaq/data",0777);
	
	srand((unsigned)time(0));
        
	cout<<"\nGenerating Public Key ...\n";	
	timer.start();// Start the timer
	cycles = rdtsc();	
	generate_matrix_MPK();
	cycles = rdtsc() - cycles;  
	timer.stop();// Stop the timer
	cout<<"Generated\n";
	cout<<"Time Elapsed in Seconds for Generating Public Key: "<<timer.getElapsedTimeInSec(); 
	cout<<"\nNumber of Cpu Cycles used for Generating Public Key: "<<cycles;
	
	cout<<"\n\nWriting Matrix B,Matrix MT,Matrix MF and Public Key in matrix(MPK) and polynomial form in \nhome\\ishtiaq\\cyclicuov_output\\Public Key.txt...\n";
	print_publickey();//It writes data in text file
	cout<<"Completed\n";
	
	cout<<"\nWriting Private Key in polynomial form in \nhome\\ishtiaq\\cyclicuov_output\\Private Key.txt...\n";
	print_privatekey();
	cout<<"Completed\n";
	
	store_privatekey();
	store_publickey();
       
	
y:	return 0;
}




 


