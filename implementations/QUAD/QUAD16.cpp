#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "Timer1.h"
using namespace std;

int choice;
int klength;
int nvar;
int L;
int D;
int*** P;
int*** Q;
int** PK;
int** QK;
int* state;
int* keystream;
int* key;
int* Gamma;
int* Delta;
int* mon;

int lookup[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                  10, 11, 12, 13, 14, 15, 0, 2, 4, 6, 8, 10, 12, 14, 3, 1, 7, 5, 11, 9, 15, 13, 0,
                  3, 6, 5, 12, 15, 10, 9, 11, 8, 13, 14, 7, 4, 1, 2, 0, 4, 8, 12, 3, 7, 11, 15, 6,
                  2, 14, 10, 5, 1, 13, 9, 0, 5, 10, 15, 7, 2, 13, 8, 14, 11, 4, 1, 9, 12, 3, 6, 0,
                  6, 12, 10, 11, 13, 7, 1, 5, 3, 9, 15, 14, 8, 2, 4, 0, 7, 14, 9, 15, 8, 1, 6, 13,
                  10, 3, 4, 2, 5, 12, 11, 0, 8, 3, 11, 6, 14, 5, 13, 12, 4, 15, 7, 10, 2, 9, 1, 0,
                  9, 1, 8, 2, 11, 3, 10, 4, 13, 5, 12, 6, 15, 7, 14, 0, 10, 7, 13, 14, 4, 9, 3,
                  15, 5, 8, 2, 1, 11, 6, 12, 0, 11, 5, 14, 10, 1, 15, 4, 7, 12, 2, 9, 13, 6, 8, 3,
                  0, 12, 11, 7, 5, 9, 14, 2, 10, 6, 1, 13, 15, 3, 4, 8, 0, 13, 9, 4, 1, 12, 8, 5,
                  2, 15, 11, 6, 3, 14, 10, 7, 0, 14, 15, 1, 13, 3, 2, 12, 9, 7, 6, 8, 4, 10, 11,
                  5, 0, 15, 13, 2, 9, 6, 4, 11, 1, 14, 12, 3, 8, 7, 5, 10 };


unsigned long long rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}
double* ave_max_min(double time[],int numberofruns)
{


	static double result_time[3];
	result_time[0]=0;
	result_time[1]=0;
	result_time[2]=10000000;


	for(int i=0;i<numberofruns;i++)
        {
	  result_time[0]=result_time[0] + time[i];

	  if(result_time[1]<=time[i])
		result_time[1]=time[i];


	  if(result_time[2]>=time[i])
		result_time[2]=time[i];

	}
	result_time[0]=result_time[0]/(double)numberofruns;
	return result_time;


}
unsigned long long* ave_max_min_cpu(unsigned long long cpucycles[],int numberofruns)
{
	 static unsigned long long result[3];
	 result[0]=0;
	 result[1]=0;
	 result[2]=1000000000;

	for(int i=0;i<numberofruns;i++)
        {
	  result[0]=result[0] + cpucycles[i];

	  if(result[1]<=cpucycles[i])
		result[1]=cpucycles[i];


	  if(result[2]>=cpucycles[i])
		result[2]=cpucycles[i];

	}
	result[0]=result[0]/numberofruns;

	return result;

}

void generate_MPK()
{
	PK= new int *[nvar];
	for (int i=0; i<nvar; i++)
		PK[i]=new int [D];


	for (int i=0; i< nvar; i++){
		for (int j=0; j<D; j++)
			PK[i][j]=(rand()%16);
	}

	QK= new int *[nvar];
		for (int i=0; i<nvar; i++)
			QK[i]=new int [D];


		for (int i=0; i< nvar; i++){
			for (int j=0; j<D; j++)
				QK[i][j]=(rand()%16);
		}

}
void generate_MP()
{

P = new int **[nvar];
for (int i = 0; i < nvar; i++) {
    		P[i] = new int*[nvar+1];
		for (int j = 0; j < (nvar+1); j++)
      			P[i][j] = new int[nvar+1];

	}
for (int i=0;i<nvar;i++){
	for (int j=0;j<nvar; j++){
		for (int k=j; k< nvar; k++){
			P[i][j][k]=(rand()%16);
		}
	}
	for (int j=0;j<nvar; j++){
			for (int k=0; k< j; k++){
				P[i][j][k]=0;
			}
		}
}

Q = new int **[nvar];
for (int i = 0; i < nvar; i++) {
    		Q[i] = new int*[nvar+1];
		for (int j = 0; j < (nvar+1); j++)
      			Q[i][j] = new int[nvar+1];

	}
for (int i=0;i<nvar;i++){
	for (int j=0;j<nvar; j++){
		for (int k=j; k< nvar; k++){
			Q[i][j][k]=(rand()%16);
		}
	}
	for (int j=0;j<nvar; j++){
			for (int k=0; k< j; k++){
				Q[i][j][k]=0;
			}
		}
}
}


void generate_Gamma()
{
	Gamma= new int[nvar];
	Delta= new int[nvar];

	for (int i=0; i< nvar;i++){
		Gamma[i]=(rand()%16);
		Delta[i]=(rand()%16);
	}

}


void generate_MPLRS()
{
P = new int **[nvar];
for (int i = 0; i < nvar; i++) {
    		P[i] = new int*[nvar+1];
		for (int j = 0; j < (nvar+1); j++)
      			P[i][j] = new int[nvar+1];

	}
for (int k=0;k<nvar;k++){
	P[k][0][0]=1;
	for (int j=1;j<nvar; j++)
		P[k][0][j]= lookup[Gamma[k]*16+P[k][0][j-1]];
	for (int i=1; i< nvar; i++){
			P[k][i][i]= lookup[Gamma[k] *16 + P[k][i-1][nvar-1]];
		for (int j =i+1; j<nvar; j++)
			P[k][i][j]=lookup[Gamma[k]*16+P[k][i][j-1]];
	}
	for (int i=0;i<nvar; i++){
			for (int j=0; j< i; j++){
				P[k][i][j]=0;
			}
		}
}

Q = new int **[nvar];
for (int i = 0; i < nvar; i++) {
    		Q[i] = new int*[nvar+1];
		for (int j = 0; j < (nvar+1); j++)
      			Q[i][j] = new int[nvar+1];

	}
for (int k=0;k<nvar;k++){
	Q[k][0][0]=1;
	for (int j=1;j<nvar; j++)
		Q[k][0][j]= lookup[Delta[k]*16+Q[k][0][j-1]];
	for (int i=1; i< nvar; i++){
			Q[k][i][i]= lookup[Delta[k] *16 + Q[k][i-1][nvar-1]];
		for (int j =i+1; j<nvar; j++)
			Q[k][i][j]=lookup[Delta[k]*16+Q[k][i][j-1]];
	}
	for (int i=0;i<nvar; i++){
			for (int j=0; j< i; j++){
				Q[k][i][j]=0;
			}
		}
}
}

void generate_MPcyclic()
{

P = new int **[nvar];
	for (int i = 0; i < nvar; i++) {
    		P[i] = new int*[nvar+1];
		for (int j = 0; j < (nvar+1); j++)
      			P[i][j] = new int[nvar+1];
	}

//first polynomial
	for (int i=0;i<nvar; i++){
		for (int j=i; j< nvar; j++){
			P[0][i][j]=(rand()%16);
		}
	}

	for (int i=0;i<nvar; i++){
		for (int j=0; j< i; j++)
				P[0][i][j]=0;
	}


//polynomials 2 - n
	for (int k=1; k<nvar; k++){
		P[k][0][0]= P[k-1][nvar-1][nvar-1];
		for (int j=1; j<nvar; j++)
			P[k][0][j]=P[k-1][0][j-1];
		for (int i=1; i<nvar; i++) {
			P[k][i][i]=P[k-1][i-1][nvar-1];
			for (int j=i+1; j<nvar; j++)
				P[k][i][j]=P[k-1][i][j-1];
		}
		for (int i=0;i<nvar; i++){
			for (int j=0; j< i; j++){
				P[k][i][j]=0;
			}
		}
	}

	Q = new int **[nvar];
	for (int i = 0; i < nvar; i++) {
    	Q[i] = new int*[nvar+1];
			for (int j = 0; j < (nvar+1); j++)
      			Q[i][j] = new int[nvar+1];
	}

	// first polynomial
	Q[0][0][0]=P[nvar-1][nvar-1][nvar-1];
	for (int j=1; j<nvar; j++)
		Q[0][0][j]=P[nvar-1][0][j-1];
	for (int i=1;i<nvar; i++){
		Q[0][i][i]=P[nvar-1][i-1][nvar-1];
		for (int j=i+1; j< nvar; j++){
			Q[0][i][j]=P[nvar-1][i][j-1];
		}
	}

	for (int i=0;i<nvar; i++){
		for (int j=0; j< i; j++)
			Q[0][i][j]=0;
	}

	//polynomials 2 - n
	for (int k=1; k<nvar; k++){
		Q[k][0][0]=Q[k-1][nvar-1][nvar-1];
		for (int j=1; j<nvar; j++)
			Q[k][0][j]=Q[k-1][0][j-1];
		for (int i=1; i<nvar; i++) {
			Q[k][i][i]=Q[k-1][i-1][nvar-1];
			for (int j=i+1; j<nvar; j++)
				Q[k][i][j]=Q[k-1][i][j-1];
		}
		for (int i=0;i<nvar; i++){
			for (int j=0; j< i; j++)
				Q[k][i][j]=0;
		}
	}
}

void print_MP(int*** MP)
{
for (int k=0; k< nvar; k++){
	for (int i=0; i< nvar; i++){
		for (int j=0; j<nvar; j++){
			cout << MP[k][i][j] << "\t";
		}
	cout << endl;}
	cout << endl;}
}

void print_MPK(int** MPK)
{
for (int k=0; k< nvar; k++){
	for (int i=0; i< D; i++){
			cout << MPK[k][i] << "\t";
		}
	cout << endl;
}
}

int* evaluatePK(int** MPK, int* state){

	int* res = new int [nvar];
	int* mon = new int [D];
	int counter =0;

	for (int j=0; j<nvar; j++){ // mon
				for (int l=j; l<nvar; l++){
					mon[counter]=lookup[state[j]*16+state[l]];
					counter++;
				}
			}

	for (int i=0; i<nvar; i++){
		res[i]=0;
		for (int j=0; j<D; j++)
			res[i]=res[i] ^ lookup[MPK[i][j]*16+mon[j]];
	}

	return res;
}

int* evaluate(int*** MP, int* state )
{
	int* temp= new int [nvar];
	int* res = new int [nvar];

for (int i=0; i< nvar; i++){
  	res[i]=0;
	for (int j=0; j<nvar; j++){
		temp[j]=0;
		for (int k=0; k<=j; k++){
			temp[j]=temp[j] ^ lookup[MP[i][k][j]*16+state[k]];
		}
	}
	for (int j=0; j<nvar; j++){
		res[i]=res[i]^lookup[temp[j]*16+state[j]];
		}
}
return res;
}

int* evaluateLRS(int*** MP, int* Beta, int* state )
{
	int* temp= new int [nvar];
	int* res = new int [nvar];

for (int i=0; i< nvar; i++){
  	res[i]=0;
  	temp[0]=state[0];
	for (int j=1; j<nvar; j++){
		temp[j]=lookup[Beta[i]*16+temp[j-1]] ^ lookup[MP[i][j][j]*16+state[j]];
	}
	for (int j=0; j<nvar; j++){
		res[i]=res[i]^lookup[temp[j]*16+state[j]];
		}
}
return res;
}

int* evaluatecyclic(int*** MP, int*** MQ, int* state )
{
	int* temp= new int [nvar];
	int* res = new int [2*nvar];

// first polynomial
	res[0]=0;
	for (int i=0; i<nvar; i++){
		temp[i]=0;
		for (int j=0; j<=i; j++){
			temp[i]=temp[i] ^ lookup[MP[0][j][i]*16+state[j]];
		}
	}

	for (int i=0; i<nvar; i++){
		res[0]=res[0] ^ lookup[temp[i]*16+state[i]];
	}

//polynomials 2 - n
for (int k=1; k< nvar; k++){
	res[k]=0;
	for (int i=(nvar-1); i>0; i--){
		temp[i]=temp[i-1] ^ lookup[MP[k][i][i]*16+state[i]];
	}
	temp[0]= lookup[MP[k][0][0]*16+state[0]];
	for (int i=0; i<nvar; i++){
		res[k]=res[k]^lookup[temp[i]*16+state[i]];
	}
}

//Q
for (int k=nvar; k< 2*nvar; k++){
	res[k]=0;
	for (int i=(nvar-1); i>0; i--){
		temp[i]=temp[i-1] ^ lookup[MQ[k-nvar][i][i]*16+state[i]];
	}
	temp[0]= lookup[MQ[k-nvar][0][0]*16+state[0]];
	for (int i=0; i<nvar; i++){
		res[k]=res[k]^lookup[temp[i]*16+state[i]];
	}
}

return res;
}

void gen_keystream()
{
	keystream= new int[nvar];

for (int i=0; i<klength; i++)
{
	keystream = evaluate(Q,state);
	/*for (int k=0; k<nvar; k++) // print keystream
		cout << keystream[k] << ", " ;
	cout << endl;*/
	state=evaluate(P,state);
}
}

void gen_keystreamPK()
{
	keystream= new int[nvar];
	mon = new int [D];


for (int i=0; i<klength; i++)
{

	keystream = evaluatePK(QK,state);
	/*for (int k=0; k<nvar; k++) // print keystream
			cout << keystream[k] << ", " ;
		cout << endl;*/
	state = evaluatePK(PK,state);
}
}

void gen_keystreamLRS()
{
	keystream= new int[nvar];

for (int i=0; i<klength; i++)
{
	keystream = evaluateLRS(Q,Delta,state);
	/*for (int k=0; k<nvar; k++) // print keystream
		cout << keystream[k] << ", " ;
	cout << endl;*/
	state=evaluateLRS(P,Gamma,state);
}
}

void gen_keystreamcyclic()
{
	int* longres= new int[2*nvar];
	keystream= new int[nvar];

for (int i=0; i<klength; i++)
{
	longres=evaluatecyclic(P,Q,state);

	/*for (int k=0; k<nvar; k++) // print keystream
		cout << longres[k+nvar] << ", " ;
	cout << endl;*/
	for (int k=0; k< nvar; k++)
		state[k]= longres[k];
}
}

int main() {
	Timer1 timer;
	srand((unsigned)time(0));
	const int numberofruns=100;
	double time_measure[numberofruns];
	div_t divresult;
		double *result_time=0;
	        unsigned long long cpu_cycles[numberofruns];
	        unsigned long long* result_cpu;
	       // unsigned long long cycles;
	cout << "QUAD over GF(16)" << endl;
	cout << "Enter 1 for QUAD using the standard approach," << endl << "2 for QUAD using the alternative approach," << endl << "3 for cyclicQUAD and 4 for QUADLRS" << endl;
	cin >> choice;
	cout << "number of variables:";
	cin >> nvar;
	cout << "Length of keystream in byte:";
	cin >> L;
	divresult= div (2*L,nvar);
	klength= divresult.quot;
	if (divresult.rem > 0)
	klength=klength+1;
	cout << numberofruns << " repetitions" << endl;
	cout << "klength: " << klength << endl;
	D= (nvar+1)*nvar/2;
	state = new int[nvar];
	for (int i=0; i<nvar; i++)
		state[i]=(rand()%16);
	if(choice==1){
	cout << endl << "***QUAD_standard approach***" << endl;
	generate_MPK();
	}
if(choice==2)
	{cout << endl << "***QUAD_alternative approach***" << endl;
	generate_MP();}
if (choice==3)
	{cout << endl << "***cyclicQUAD***" << endl;
	generate_MPcyclic();}
if (choice==4){
	cout << endl << "***QUADLRS***" << endl;
	generate_Gamma();
	generate_MPLRS();
	}

if (choice==1){
	for(int i=0;i<numberofruns;i++)
		{
		timer.start();// Start the timer
		cpu_cycles[i] = rdtsc();
		gen_keystreamPK();
		//print_keystream();
	    cpu_cycles[i] = rdtsc() - cpu_cycles[i];
		timer.stop();// Stop the timer
	    time_measure[i]=timer.getElapsedTimeInMicroSec();
	    }
}

if (choice==2){
	for(int i=0;i<numberofruns;i++)
		{
		timer.start();// Start the timer
		cpu_cycles[i] = rdtsc();
		gen_keystream();
		//print_keystream();
	    cpu_cycles[i] = rdtsc() - cpu_cycles[i];
		timer.stop();// Stop the timer
	    time_measure[i]=timer.getElapsedTimeInMicroSec();
	    }
}

if (choice==3){
	for(int i=0;i<numberofruns;i++)
		{
		timer.start();// Start the timer
		cpu_cycles[i] = rdtsc();
		gen_keystreamcyclic();
		//print_keystream();
	    cpu_cycles[i] = rdtsc() - cpu_cycles[i];
		timer.stop();// Stop the timer
	    time_measure[i]=timer.getElapsedTimeInMicroSec();
	    }
}

if (choice==4){
	for(int i=0;i<numberofruns;i++)
		{
		timer.start();// Start the timer
		cpu_cycles[i] = rdtsc();
		gen_keystreamLRS();
		//print_keystream();
	    cpu_cycles[i] = rdtsc() - cpu_cycles[i];
		timer.stop();// Stop the timer
	    time_measure[i]=timer.getElapsedTimeInMicroSec();
		}
}
		result_time=ave_max_min(time_measure,numberofruns);
		cout<<"Time Elapsed in  Micro-Seconds for key stream generation\n";
		cout<<"Average Time: "<<result_time[0]<<"\nMaximum Time: "<<result_time[1]<<"\nMinimum Time: "<<result_time[2];

		result_cpu=ave_max_min_cpu(cpu_cycles,numberofruns);
		cout<<"\n\nNumber of Cpu Cycles used for key stream generation\n";
		cout<<"Average Cpu Cycles: "<<result_cpu[0]<<"\nMaximum Cpu Cycles: "<<result_cpu[1]<<"\nMinimum Cpu Cycles: "<<result_cpu[2];

		return 0;}
