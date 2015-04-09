// Ishtiaqcourage.cpp : Defines the entry point for the console application.
//

#include "basic.h"
#include "sign.h"
#include "loadprivatekey.h"
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
        
	
	cout<<"\n**********Implementation of Signature Generation of Rainbow Scheme***********\n";


	
	srand((unsigned)time(0));
	
    	load_privatekey();
	
	cout<<"\nGenerating Signature ...\n";

	for(int i=0;i<length;i++)
	{   
		timer.start();// Start the timer
		cpu_cycles[i] = rdtsc(); 
    		generatesignature();
                cpu_cycles[i] = rdtsc() - cpu_cycles[i];         
		timer.stop();// Stop the timer
                time_measure[i]=timer.getElapsedTimeInMicroSec();
		
	}

	cout<<"Generated\n";

	result_time=ave_max_min(time_measure,length);
	cout<<"Time Elapsed in  Micro-Seconds for Generating Signature\n";
	cout<<"Average Time: "<<result_time[0]<<"\nMaximum Time: "<<result_time[1]<<"\nMinimum Time: "<<result_time[2];
	result_cpu=ave_max_min_cpu(cpu_cycles,length);
	cout<<"\n\nNumber of Cpu Cycles used for Generating Signature\n";
	cout<<"Average Cpu Cycles: "<<result_cpu[0]<<"\nMaximum Cpu Cycles: "<<result_cpu[1]<<"\nMinimum Cpu Cycles: "<<result_cpu[2];	

	cout<<"\nWriting the following Hashvalues, Vinegar variables, Oil Variables(O1 and O2)\nand Signature in home\\ishtiaq\\rainbow_output\\Signature Generation.txt...\n";
	print_sign_generation();
	cout<<"Completed\n";

	store_hash_sign();
	
	
	return 0;
}




 


