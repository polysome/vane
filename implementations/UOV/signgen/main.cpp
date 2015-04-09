// Ishtiaqcourage.cpp : Defines the entry point for the console application.
//

#include "basic.h"




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
        
	
    cout<<"\n**********Implementation of Unbalanced Oil and Vinegar Scheme***********\n";
	


	//n=oilvariables+vinegarvariables;
	//D = (vinegarvariables*(2*oilvariables+vinegarvariables+1))/2;
	
    mkdir("/home/ishtiaq/uov_output",0777);
	
	srand((unsigned)time(0));
        
	load_privatekey();
	
	cout<<"\n\nGenerating Signature ...\n";
	
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
	
	
	
	
	
	
        	
    cout<<"\nWriting hashvalues,oil and Vinegar Variables and Signature in \nhome\\ishtiaq\\uov_output\\Signature Generation.txt...\n";
	print_sign_generation();
	cout<<"Completed\n";
	store_hash_sign();
	
	
	
	
	
	return 0;
}




 





